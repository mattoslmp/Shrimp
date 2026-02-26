from __future__ import annotations

import csv
import gzip
import json
import os
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import requests
import streamlit as st
from Bio import Entrez, SeqIO
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


# ==============================
# Configurações gerais
# ==============================
APP_TITLE = "Shrimp 16S Explorer"
DATA_DIR = Path("data")
DOWNLOAD_DIR = DATA_DIR / "downloads"
REPORTS_DIR = DATA_DIR / "reports"
QIIME_DIR = DATA_DIR / "qiime2"
for d in [DATA_DIR, DOWNLOAD_DIR, REPORTS_DIR, QIIME_DIR]:
    d.mkdir(parents=True, exist_ok=True)


@dataclass
class SearchConfig:
    species: str = "Penaeus vannamei"
    marker: str = "16S"
    max_records: int = 200


class HTTPClient:
    """Sessão HTTP com retry para ENA/NCBI APIs."""

    @staticmethod
    def build_session() -> requests.Session:
        session = requests.Session()
        retry = Retry(
            total=5,
            backoff_factor=1,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["GET", "POST"],
        )
        session.mount("https://", HTTPAdapter(max_retries=retry))
        session.mount("http://", HTTPAdapter(max_retries=retry))
        return session


class NCBIService:
    """Consulta registros NCBI/SRA relacionados a 16S e espécie alvo."""

    def __init__(self, email: str, api_key: str = "") -> None:
        Entrez.email = email.strip() or "anonymous@example.org"
        if api_key.strip():
            Entrez.api_key = api_key.strip()

    def search_sra_ids(self, config: SearchConfig) -> List[str]:
        term = (
            f'("{config.species}"[Organism]) AND ({config.marker}[All Fields]) '
            "AND (biomol genomic[Properties] OR amplicon[All Fields])"
        )
        with Entrez.esearch(db="sra", term=term, retmax=config.max_records) as handle:
            record = Entrez.read(handle)
        return record.get("IdList", [])

    def sra_summaries(self, id_list: List[str]) -> pd.DataFrame:
        if not id_list:
            return pd.DataFrame()
        chunks = [id_list[i : i + 200] for i in range(0, len(id_list), 200)]
        rows: List[Dict[str, str]] = []

        for chunk in chunks:
            with Entrez.esummary(db="sra", id=",".join(chunk), retmode="json") as handle:
                summary = json.load(handle)
            result = summary.get("result", {})
            uids = result.get("uids", [])
            for uid in uids:
                item = result.get(uid, {})
                exp_xml = item.get("expxml", "")
                runs = item.get("runs", "")
                rows.append(
                    {
                        "uid": uid,
                        "title": item.get("title", ""),
                        "organism": item.get("organism", ""),
                        "bioproject": item.get("bioproject", ""),
                        "biosample": item.get("biosample", ""),
                        "experiment_xml": exp_xml,
                        "runs_raw": runs,
                    }
                )
            time.sleep(0.34)

        df = pd.DataFrame(rows)
        df["run_accession"] = df["runs_raw"].str.extract(r"(SRR\d+|ERR\d+|DRR\d+)")
        return df


class ENAService:
    """Obtém links de FASTQ na ENA via accession de run."""

    BASE_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"

    def __init__(self) -> None:
        self.session = HTTPClient.build_session()

    def fetch_fastq_links(self, run_accession: str) -> Dict[str, str]:
        params = {
            "accession": run_accession,
            "result": "read_run",
            "fields": "run_accession,fastq_ftp,fastq_md5,library_layout,instrument_platform,scientific_name",
            "format": "json",
        }
        r = self.session.get(self.BASE_URL, params=params, timeout=60)
        r.raise_for_status()
        payload = r.json()
        if not payload:
            return {"run_accession": run_accession, "fastq_ftp": "", "fastq_md5": ""}
        row = payload[0]
        return {
            "run_accession": row.get("run_accession", run_accession),
            "fastq_ftp": row.get("fastq_ftp", ""),
            "fastq_md5": row.get("fastq_md5", ""),
            "library_layout": row.get("library_layout", ""),
            "instrument_platform": row.get("instrument_platform", ""),
            "scientific_name": row.get("scientific_name", ""),
        }


class Downloader:
    """Baixa FASTQ por URL http/https/ftp."""

    def __init__(self) -> None:
        self.session = HTTPClient.build_session()

    def normalize_ena_link(self, link: str) -> str:
        if link.startswith("ftp.sra.ebi.ac.uk"):
            return f"https://{link}"
        if link.startswith("ftp://"):
            return link.replace("ftp://", "https://")
        return link

    def download_file(self, url: str, out_dir: Path) -> Path:
        out_dir.mkdir(parents=True, exist_ok=True)
        fixed = self.normalize_ena_link(url)
        name = fixed.rstrip("/").split("/")[-1]
        out_path = out_dir / name

        with self.session.get(fixed, stream=True, timeout=300) as r:
            r.raise_for_status()
            with open(out_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        f.write(chunk)
        return out_path


class QualityAnalyzer:
    """Métricas básicas de qualidade para FASTQ."""

    def analyze_fastq(self, file_path: Path, max_reads: int = 50000) -> Dict[str, float]:
        opener = gzip.open if file_path.suffix == ".gz" else open
        n_reads = 0
        lengths: List[int] = []
        mean_q: List[float] = []

        with opener(file_path, "rt", encoding="utf-8", errors="ignore") as handle:
            for rec in SeqIO.parse(handle, "fastq"):
                q = rec.letter_annotations.get("phred_quality", [])
                if not q:
                    continue
                n_reads += 1
                lengths.append(len(rec.seq))
                mean_q.append(sum(q) / len(q))
                if n_reads >= max_reads:
                    break

        if n_reads == 0:
            return {
                "file": file_path.name,
                "reads_sampled": 0,
                "mean_length": 0,
                "mean_quality": 0,
                "q30_rate": 0,
            }

        q30 = sum(1 for v in mean_q if v >= 30) / len(mean_q)
        return {
            "file": file_path.name,
            "reads_sampled": n_reads,
            "mean_length": round(sum(lengths) / len(lengths), 2),
            "mean_quality": round(sum(mean_q) / len(mean_q), 2),
            "q30_rate": round(q30, 4),
        }


class QIIME2Pipeline:
    """Gera manifest e executa pipeline QIIME2 (opcional)."""

    def __init__(self, qiime_env: str = "qiime2") -> None:
        self.qiime_env = qiime_env

    def build_manifest(self, fastq_files: List[Path], manifest_path: Path) -> Path:
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        rows: List[Tuple[str, str, str]] = []

        for f in fastq_files:
            sample = f.name.split("_")[0]
            lower = f.name.lower()
            if "_r1" in lower or "_1" in lower:
                direction = "forward"
            elif "_r2" in lower or "_2" in lower:
                direction = "reverse"
            else:
                direction = "forward"
            rows.append((sample, str(f.resolve()), direction))

        with open(manifest_path, "w", newline="", encoding="utf-8") as out:
            writer = csv.writer(out)
            writer.writerow(["sample-id", "absolute-filepath", "direction"])
            writer.writerows(rows)

        return manifest_path

    def commands(self, manifest_path: Path, out_dir: Path) -> List[str]:
        out_dir.mkdir(parents=True, exist_ok=True)
        return [
            f"qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path {manifest_path} --output-path {out_dir / 'demux.qza'} --input-format PairedEndFastqManifestPhred33V2",
            f"qiime demux summarize --i-data {out_dir / 'demux.qza'} --o-visualization {out_dir / 'demux.qzv'}",
            f"qiime dada2 denoise-paired --i-demultiplexed-seqs {out_dir / 'demux.qza'} --p-trunc-len-f 240 --p-trunc-len-r 200 --o-table {out_dir / 'table.qza'} --o-representative-sequences {out_dir / 'rep-seqs.qza'} --o-denoising-stats {out_dir / 'stats.qza'}",
            f"qiime feature-table summarize --i-table {out_dir / 'table.qza'} --o-visualization {out_dir / 'table.qzv'}",
            f"qiime feature-table tabulate-seqs --i-data {out_dir / 'rep-seqs.qza'} --o-visualization {out_dir / 'rep-seqs.qzv'}",
        ]

    def run(self, commands: List[str]) -> List[Tuple[str, int, str]]:
        logs: List[Tuple[str, int, str]] = []
        for cmd in commands:
            proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            logs.append((cmd, proc.returncode, proc.stdout + "\n" + proc.stderr))
            if proc.returncode != 0:
                break
        return logs


# ==============================
# UI
# ==============================

def apply_theme() -> None:
    st.set_page_config(page_title=APP_TITLE, layout="wide", page_icon="🦐")
    st.markdown(
        """
        <style>
        .stApp {background: linear-gradient(180deg, #041c32 0%, #064663 45%, #0b8f8c 100%);} 
        [data-testid="stHeader"] {background: transparent;}
        .block-container {padding-top: 1.2rem;}
        .hero-card {
            background: rgba(255,255,255,0.08);
            border: 1px solid rgba(255,255,255,0.18);
            border-radius: 18px;
            padding: 18px;
            backdrop-filter: blur(8px);
            color: #e8f7ff;
        }
        .metric-box {
            background: rgba(255,255,255,0.92);
            border-radius: 12px;
            padding: 8px 14px;
            border-left: 6px solid #0b8f8c;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )


def enrich_with_ena_links(df: pd.DataFrame, ena: ENAService) -> pd.DataFrame:
    if df.empty or "run_accession" not in df.columns:
        return df

    rows = []
    for run in df["run_accession"].dropna().unique().tolist():
        try:
            rows.append(ena.fetch_fastq_links(run))
        except Exception as e:
            rows.append({"run_accession": run, "fastq_ftp": "", "error": str(e)})

    links_df = pd.DataFrame(rows)
    return df.merge(links_df, on="run_accession", how="left")


def user_mode() -> None:
    st.markdown("<div class='hero-card'><h2>🧬 Coleta de dados 16S para camarão</h2><p>Busca automática no NCBI + ENA com download de FASTQ e metadados.</p></div>", unsafe_allow_html=True)

    email = st.text_input("Email NCBI", value=os.getenv("NCBI_EMAIL", ""))
    api_key = st.text_input("API Key NCBI (opcional)", value=os.getenv("NCBI_API_KEY", ""), type="password")
    species = st.text_input("Espécie alvo", value="Penaeus vannamei")
    marker = st.text_input("Marcador", value="16S")
    max_records = st.slider("Máximo de registros SRA", min_value=20, max_value=1000, value=200, step=20)

    if st.button("🔎 Buscar datasets"):
        try:
            ncbi = NCBIService(email=email, api_key=api_key)
            ena = ENAService()
            cfg = SearchConfig(species=species, marker=marker, max_records=max_records)

            ids = ncbi.search_sra_ids(cfg)
            base_df = ncbi.sra_summaries(ids)
            data = enrich_with_ena_links(base_df, ena)
            st.session_state["datasets"] = data

            out_csv = REPORTS_DIR / f"datasets_{species.replace(' ', '_')}.csv"
            data.to_csv(out_csv, index=False)

            c1, c2, c3 = st.columns(3)
            c1.metric("SRA IDs", len(ids))
            c2.metric("Runs detectadas", data["run_accession"].nunique() if not data.empty else 0)
            c3.metric("Com link FASTQ", int(data["fastq_ftp"].fillna("").str.len().gt(0).sum()) if not data.empty else 0)

            st.success(f"Metadados salvos em: {out_csv}")
            st.dataframe(data, use_container_width=True)
        except Exception as e:
            st.error(f"Falha na busca: {e}")

    if "datasets" in st.session_state and not st.session_state["datasets"].empty:
        df = st.session_state["datasets"]
        st.subheader("Download de FASTQ")
        only_with_link = df[df["fastq_ftp"].fillna("") != ""].copy()
        st.write(f"Runs com links disponíveis: {len(only_with_link)}")

        if st.button("⬇️ Baixar primeiros 10 arquivos FASTQ"):
            dl = Downloader()
            paths = []
            for _, row in only_with_link.head(10).iterrows():
                links = str(row.get("fastq_ftp", "")).split(";")
                for link in links:
                    link = link.strip()
                    if not link:
                        continue
                    try:
                        p = dl.download_file(link, DOWNLOAD_DIR)
                        paths.append(p)
                    except Exception as e:
                        st.warning(f"Erro no download {link}: {e}")
            st.session_state["downloaded"] = paths
            st.success(f"Arquivos baixados: {len(paths)}")

    if "downloaded" in st.session_state and st.session_state["downloaded"]:
        st.subheader("Quality check básico")
        qa = QualityAnalyzer()
        rows = []
        for fp in st.session_state["downloaded"]:
            try:
                rows.append(qa.analyze_fastq(fp))
            except Exception as e:
                rows.append({"file": fp.name, "reads_sampled": 0, "mean_length": 0, "mean_quality": 0, "q30_rate": 0, "error": str(e)})
        qdf = pd.DataFrame(rows)
        st.dataframe(qdf, use_container_width=True)
        qfile = REPORTS_DIR / "quality_report.csv"
        qdf.to_csv(qfile, index=False)
        st.info(f"Relatório de qualidade salvo em {qfile}")


def admin_mode() -> None:
    st.markdown("<div class='hero-card'><h2>🛠️ Modo Admin</h2><p>Pipeline QIIME2, manifest e execução supervisionada.</p></div>", unsafe_allow_html=True)

    pwd = st.text_input("Senha admin", type="password")
    expected = os.getenv("ADMIN_PASSWORD", "admin123")
    if pwd != expected:
        st.warning("Informe a senha admin para liberar ações avançadas.")
        return

    st.success("Admin autenticado.")
    files = sorted(DOWNLOAD_DIR.glob("*.fastq*"))
    st.write(f"FASTQ detectados: {len(files)}")

    if files:
        st.code("\n".join(str(f) for f in files[:30]))

    qiime_env = st.text_input("Conda env QIIME2", value="qiime2")
    run_pipeline = st.checkbox("Executar comandos QIIME2 agora", value=False)

    if st.button("⚙️ Gerar manifest e comandos"):
        pipeline = QIIME2Pipeline(qiime_env=qiime_env)
        manifest = pipeline.build_manifest(files, QIIME_DIR / "manifest.csv")
        cmds = pipeline.commands(manifest, QIIME_DIR)

        st.write("Manifest:", manifest)
        st.code("\n".join(cmds), language="bash")

        if run_pipeline:
            logs = pipeline.run(cmds)
            for cmd, rc, log in logs:
                if rc == 0:
                    st.success(f"OK: {cmd}")
                else:
                    st.error(f"Falhou ({rc}): {cmd}")
                st.text_area("Log", value=log, height=180)


def main() -> None:
    apply_theme()
    st.title("🦐 Shrimp 16S Explorer")
    st.caption("Coleta de metadados e reads de SRA/ENA, QC e preparação para QIIME2")

    mode = st.sidebar.radio("Modo", ["User", "Admin"], index=0)
    st.sidebar.info(
        "Use o modo User para obter metadados e downloads; "
        "use Admin para gerar manifest e rodar QIIME2."
    )

    if mode == "User":
        user_mode()
    else:
        admin_mode()


if __name__ == "__main__":
    main()
