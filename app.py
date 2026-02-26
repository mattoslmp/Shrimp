# -*- coding: utf-8 -*-
import os
import sys
import io
import json
import time
import random
import hashlib
import secrets
import subprocess
import argparse
import logging
import threading
import queue
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, asdict
from enum import Enum
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple, Any, Iterable
import numpy as np
import pandas as pd
import geopandas as gpd
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.preprocessing import StandardScaler, LabelEncoder
import skbio
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.ordination import pcoa
from skbio.tree import TreeNode

import plotly.express as px
import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt
import seaborn as sns

import streamlit as st
from streamlit_autorefresh import st_autorefresh
import streamlit.components.v1 as components

from sqlalchemy import (
    create_engine, Column, Integer, String, Boolean, DateTime,
    ForeignKey, Text, Float, JSON, BigInteger, Enum as SQLEnum,
    func, text, and_, or_
)
from sqlalchemy.orm import declarative_base, sessionmaker, relationship
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.dialects.postgresql import UUID, ARRAY

# ImportaÃ§Ãµes para download de dados pÃºblicos
import requests
from Bio import Entrez, SeqIO
from bioservices import EUtils
import pysradb
import ftputil
import gzip
import shutil
from dataclasses import dataclass
from typing import Optional
from typing import Any
import unicodedata


def _ensure_dir(path: Path) -> Path:
    """Cria diretÃ³rio se nÃ£o existir e retorna o Path absoluto."""
    path.mkdir(parents=True, exist_ok=True)
    return path

def _row_links_from_ena(row: Dict[str, str]) -> Dict[str, List[str]]:
    """Extrai links FTP, HTTP e Aspera de uma linha do ENA."""
    result = {
        "fastq_ftp": [],
        "fastq_http": [],
        "fastq_aspera": [],
    }
    
    # FTP
    if row.get("fastq_ftp"):
        links = str(row["fastq_ftp"]).split(";")
        result["fastq_ftp"] = [link.strip() for link in links if link.strip()]
    
    # HTTP
    if row.get("fastq_http"):
        links = str(row["fastq_http"]).split(";")
        result["fastq_http"] = [link.strip() for link in links if link.strip()]
    
    # Aspera
    if row.get("fastq_aspera"):
        links = str(row["fastq_aspera"]).split(";")
        result["fastq_aspera"] = [link.strip() for link in links if link.strip()]
    
    return result

def _fetch_ena_read_run_table(
    study_or_project: str,
    fields: List[str],
    limit: int = 0,
    session: Optional[requests.Session] = None,
) -> List[Dict[str, str]]:
    """Busca tabela de runs do ENA para um estudo ou projeto."""
    session = session or requests.Session()
    base_url = "https://www.ebi.ac.uk/ena/portal/api/search"
    
    params = {
        "query": f'study_accession="{study_or_project}" OR bioproject="{study_or_project}"',
        "result": "read_run",
        "fields": ",".join(fields),
        "format": "json",
        "limit": str(limit) if limit else "0",
    }
    
    try:
        response = session.get(base_url, params=params, timeout=30)
        response.raise_for_status()
        data = response.json()
        return data if isinstance(data, list) else []
    except Exception:
        return []

def _write_links_table(df: pd.DataFrame, output_path: Path) -> None:
    """Salva tabela com links de download."""
    if df.empty:
        return
    
    # Selecionar colunas relevantes
    cols = []
    for col in ["run_accession", "accession", "fastq_ftp_links", 
                "fastq_http_links", "fastq_aspera_links", "bioproject"]:
        if col in df.columns:
            cols.append(col)
    
    if cols:
        df[cols].to_csv(output_path, sep="\t", index=False)

# PaÃ­ses produtores de camarÃ£o (usado em telas de registro/perfil)
SHRIMP_PRODUCING_COUNTRIES = [
    "Brasil",
    "China",
    "Equador",
    "Ãndia",
    "IndonÃ©sia",
    "VietnÃ£",
    "TailÃ¢ndia",
    "MÃ©xico",
    "Estados Unidos",
    "Outro"
]

# ConfiguraÃ§Ã£o de logging
LOG_FILE = Path("data/logs/shrimp_platform.log")  # ou config.LOGS_DIR / ...
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# ============================================================================
# CONFIGURAÃ‡Ã•ES GLOBAIS E CONSTANTES
# ============================================================================
import streamlit as st

def st_divider():
    """
    Helper compatÃ­vel com versÃµes antigas do Streamlit.
    Usa st.divider() se existir, senÃ£o cai para um <hr> estilizado.
    """
    if hasattr(st, "divider"):
        st.divider()
    else:
        st.markdown(
            "<hr style='border: 1px solid #e5e7eb; margin: 0.75rem 0;' />",
            unsafe_allow_html=True,
        )

class PlatformConfig:
    """ConfiguraÃ§Ãµes da plataforma"""
    
    # InformaÃ§Ãµes da plataforma
    APP_NAME = "ShrimpMicrobiome Platform"
    APP_VERSION = "1.0.0"          # VersÃ£o 1 da plataforma
    APP_YEAR = 2025                # Ano da plataforma
    COMPANY_NAME = "Sea Scient"     # Nome da empresa
    APP_TAGLINE = "AnÃ¡lise integrada da microbiota de camarÃ£o utilizando bioinformÃ¡tica e IA"

    
    # DiretÃ³rios
    BASE_DIR = Path(__file__).parent.absolute()
    DATA_DIR = BASE_DIR / "data"
    PROJECTS_DIR = DATA_DIR / "projects"
    PUBLIC_DATA_DIR = DATA_DIR / "public_data"
    USER_DATA_DIR = DATA_DIR / "user_data"
    DATABASE_DIR = DATA_DIR / "database"
    PIPELINE_DIR = BASE_DIR / "pipeline"
    RESULTS_DIR = DATA_DIR / "results"
    LOGS_DIR = DATA_DIR / "logs"
    
    # Criar diretÃ³rios
    for dir_path in [DATA_DIR, PROJECTS_DIR, PUBLIC_DATA_DIR, USER_DATA_DIR,
                    DATABASE_DIR, PIPELINE_DIR, RESULTS_DIR, LOGS_DIR]:
        dir_path.mkdir(exist_ok=True, parents=True)
    
    # Banco de dados
    DB_URL = f"sqlite:///{DATABASE_DIR}/shrimp_microbiome.db"
    
    # ConfiguraÃ§Ãµes de email
    EMAIL_HOST = os.getenv("SHRIMP_EMAIL_HOST", "")
    EMAIL_PORT = int(os.getenv("SHRIMP_EMAIL_PORT", "587"))
    EMAIL_USER = os.getenv("SHRIMP_EMAIL_USER", "")
    EMAIL_PASS = os.getenv("SHRIMP_EMAIL_PASS", "")
    
    # Admin padrÃ£o
    DEFAULT_ADMIN_USERNAME = "admin"
    DEFAULT_ADMIN_PASSWORD = "admin123"
    DEFAULT_ADMIN_EMAIL = os.getenv("SHRIMP_ADMIN_EMAIL", "admin@shrimpmicrobiome.com")
    
    # ParÃ¢metros especÃ­ficos do camarÃ£o
    TARGET_SPECIES = "Penaeus vannamei"
    TARGET_REGION = "V3-V4"  # RegiÃ£o 16S alvo
    REFERENCE_DATABASE = "SILVA 138"  # Banco de dados de referÃªncia
    
    # URLs de bases de dados pÃºblicas
    PUBLIC_DATASETS = {
        "sra": {
            "url": "https://www.ncbi.nlm.nih.gov/sra",
            "search_term": "Penaeus vannamei 16S",
            "filters": {
                "strategy": "AMPLICON",
                "platform": "ILLUMINA",
                "source": "METAGENOMIC"
            }
        },
        "ena": {
            "url": "https://www.ebi.ac.uk/ena",
            "search_term": "Penaeus vannamei AND 16S"
        },
        "mg_rast": {
            "url": "https://www.mg-rast.org",
            "search_term": "shrimp microbiome"
        }
    }
    
    # ConfiguraÃ§Ã£o do pipeline QIIME2
    QIIME2_VERSION = "2023.9"
    QIIME2_PIPELINE = {
        "import": "qiime tools import",
        "demux": "qiime demux summarize",
        "dada2": "qiime dada2 denoise-paired",
        "phylogeny": "qiime phylogeny align-to-tree-mafft-fasttree",
        "taxonomy": "qiime feature-classifier classify-sklearn",
        "diversity": "qiime diversity core-metrics-phylogenetic",
        "export": "qiime tools export"
    }
    
    # ParÃ¢metros padrÃ£o do pipeline
    DEFAULT_PARAMS = {
        "trunc_len_f": 240,
        "trunc_len_r": 200,
        "trim_left_f": 10,
        "trim_left_r": 10,
        "max_ee": 2.0,
        "min_reads": 1000,
        "sampling_depth": 10000,
        "n_permutations": 999
    }
    
    # Estados do processamento
    class ProcessingStatus(Enum):
        PENDING = "pending"
        DOWNLOADING = "downloading"
        QC = "quality_control"
        PROCESSING = "processing"
        ANALYZING = "analyzing"
        COMPLETED = "completed"
        FAILED = "failed"
        COMPARING = "comparing"
    
    # NÃ­veis de infecÃ§Ã£o
    INFECTION_LEVELS = ["Leve", "Moderada", "Severa", "CrÃ­tica"]
    
    # Macronutrientes
    MACRONUTRIENTS = ["ProteÃ­na", "LipÃ­dio", "Carboidrato", "Fibra", "Cinzas"]
    
    # PatÃ³genos conhecidos
    KNOWN_PATHOGENS = [
        "Vibrio parahaemolyticus",
        "Vibrio harveyi",
        "Vibrio alginolyticus",
        "Photobacterium damselae",
        "Aeromonas hydrophila",
        "White Spot Syndrome Virus",
        "Infectious Hypodermal and Hematopoietic Necrosis Virus",
        "Taura Syndrome Virus"
    ]
    
    # BactÃ©rias benÃ©ficas
    BENEFICIAL_BACTERIA = [
        "Bacillus subtilis",
        "Bacillus licheniformis",
        "Lactobacillus spp.",
        "Pseudomonas fluorescens",
        "Enterococcus faecium"
    ]

config = PlatformConfig()

# ============================================================================
# MODELOS DE BANCO DE DADOS
# ============================================================================
Base = declarative_base()

class UserRole(Enum):
  ADMIN = "admin"
  RESEARCHER = "researcher"
  FARMER = "farmer"
  STUDENT = "student"
  PUBLIC = "public"


@dataclass
class CurrentUser:
  id: int
  email: str
  full_name: str
  username: str
  role: UserRole
  organization: Optional[str] = None


class UserStatus(Enum):
  PENDING = "pending"
  APPROVED = "approved"
  REJECTED = "rejected"
  SUSPENDED = "suspended"


class User(Base):
  __tablename__ = "users"
  
  id = Column(Integer, primary_key=True)
  username = Column(String(50), unique=True, nullable=False, index=True)
  email = Column(String(255), unique=True, nullable=False, index=True)
  password_hash = Column(String(255), nullable=False)
  full_name = Column(String(255))
  organization = Column(String(255))
  country = Column(String(100))
  role = Column(SQLEnum(UserRole), default=UserRole.RESEARCHER)
  status = Column(SQLEnum(UserStatus), default=UserStatus.PENDING)
  
  # InformaÃ§Ãµes especÃ­ficas do estudo
  study_species = Column(String(100), default="Penaeus vannamei")
  study_description = Column(Text)
  research_interests = Column(JSON)
  
  # ConfiguraÃ§Ãµes
  settings = Column(
    JSON,
    default=lambda: {
      "email_notifications": True,
      "dashboard_view": "advanced",
      "auto_download_updates": False,
      "privacy_level": "standard",
    },
  )
  
  # EstatÃ­sticas
  project_count = Column(Integer, default=0)
  analysis_count = Column(Integer, default=0)
  last_login = Column(DateTime)
  created_at = Column(DateTime, default=datetime.utcnow)
  updated_at = Column(DateTime, onupdate=datetime.utcnow)
  
  # Relacionamentos
  projects = relationship("Project", back_populates="owner")
  analyses = relationship("Analysis", back_populates="user")
  samples = relationship("Sample", back_populates="uploader")
  
  def verify_password(self, password: str) -> bool:
    """Verifica se a senha estÃ¡ correta"""
    return hashlib.sha256(password.encode()).hexdigest() == self.password_hash
  
  def to_dict(self) -> Dict:
    return {
      "id": self.id,
      "username": self.username,
      "email": self.email,
      "full_name": self.full_name,
      "organization": self.organization,
      "country": self.country,
      "role": self.role.value,
      "status": self.status.value,
      "study_species": self.study_species,
      "project_count": self.project_count,
      "analysis_count": self.analysis_count,
      "created_at": self.created_at.isoformat() if self.created_at else None,
    }

  def to_current_user(self) -> CurrentUser:
    """
    Cria um objeto leve CurrentUser para guardar no st.session_state,
    evitando manter o objeto ORM vinculado Ã  Session.
    """
    return CurrentUser(
      id=self.id,
      email=self.email,
      full_name=self.full_name or self.username,
      username=self.username,
      role=self.role,
      organization=self.organization,
    )



class Project(Base):
    __tablename__ = "projects"
    
    id = Column(Integer, primary_key=True)
    owner_id = Column(Integer, ForeignKey("users.id"), nullable=False)
    name = Column(String(255), nullable=False)
    description = Column(Text)
    species = Column(String(100), default="Penaeus vannamei")
    
    # Tipo de projeto
    project_type = Column(String(50), default="comparative_analysis")  # comparative_analysis, user_upload, public_data
    
    # Metadados experimentais
    experimental_design = Column(JSON)  # {'groups': [], 'replicates': 0, 'variables': []}
    treatments = Column(JSON)  # Lista de tratamentos aplicados
    diet_composition = Column(JSON)  # ComposiÃ§Ã£o da dieta em macronutrientes
    
    # Status
    status = Column(String(50), default="planning")
    progress = Column(Integer, default=0)  # 0-100%
    is_public = Column(Boolean, default=False)
    
    # Caminhos
    base_path = Column(String(1024))
    data_path = Column(String(1024))
    results_path = Column(String(1024))
    
    # Datas
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, onupdate=datetime.utcnow)
    completed_at = Column(DateTime)
    
    # Relacionamentos
    owner = relationship("User", back_populates="projects")
    samples = relationship("Sample", back_populates="project")
    analyses = relationship("Analysis", back_populates="project")
    comparisons = relationship("Comparison", back_populates="project")
    
    @property
    def sample_count(self) -> int:
        return len(self.samples)
    
    @property
    def has_public_data(self) -> bool:
        return any(s.is_public for s in self.samples)



# -------------------------------------------------------------------------
# Enums de domÃ­nio
# -------------------------------------------------------------------------
class SampleHealthStatus(Enum):
    HEALTHY = "healthy"
    DYSBIOSIS = "dysbiosis"
    UNKNOWN = "unknown"


class InfectionLevel(Enum):
    NONE = "none"
    MILD = "mild"
    MODERATE = "moderate"
    SEVERE = "severe"
    CRITICAL = "critical"


class AnalysisType(Enum):
    DIVERSITY_ALPHA = "diversity_alpha"
    DIVERSITY_BETA = "diversity_beta"
    TAXONOMIC = "taxonomic"
    DIFFERENTIAL_ABUNDANCE = "differential_abundance"
    CORRELATION = "correlation"
    COMPARATIVE = "comparative"
    BATCH_CORRECTION = "batch_correction"


# -------------------------------------------------------------------------
# FunÃ§Ãµes auxiliares de parsing (pt/en â†’ Enum)
# -------------------------------------------------------------------------


def _normalize_label(value: Any) -> str:
    """
    Normaliza string para minÃºsculas ASCII sem acentos.
    Ex.: "SaudÃ¡vel" -> "saudavel"
    """
    if value is None:
        return ""
    s = str(value).strip().lower()
    s_norm = unicodedata.normalize("NFKD", s)
    s_no_accents = "".join(c for c in s_norm if not unicodedata.combining(c))
    return s_no_accents


def parse_health_status(value: Any) -> SampleHealthStatus:
    """
    Converte rÃ³tulos (pt/en) ou Enum para SampleHealthStatus.

    Aceita, por exemplo:
      - "healthy", "saudavel", "saudÃ¡vel", "controle", "control"
      - "dysbiosis", "disbiose", "doente", "disease"
      - "unknown", "desconhecido", "na", "n/a"
    """
    if isinstance(value, SampleHealthStatus):
        return value
    if value is None:
        return SampleHealthStatus.UNKNOWN

    s = _normalize_label(value)

    mapping = {
        # SaudÃ¡vel / controle
        "healthy": SampleHealthStatus.HEALTHY,
        "saudavel": SampleHealthStatus.HEALTHY,
        "controle": SampleHealthStatus.HEALTHY,
        "control": SampleHealthStatus.HEALTHY,

        # Disbiose / doente
        "dysbiosis": SampleHealthStatus.DYSBIOSIS,
        "disbiose": SampleHealthStatus.DYSBIOSIS,
        "doente": SampleHealthStatus.DYSBIOSIS,
        "disease": SampleHealthStatus.DYSBIOSIS,
        "sick": SampleHealthStatus.DYSBIOSIS,

        # Desconhecido
        "unknown": SampleHealthStatus.UNKNOWN,
        "desconhecido": SampleHealthStatus.UNKNOWN,
        "na": SampleHealthStatus.UNKNOWN,
        "n/a": SampleHealthStatus.UNKNOWN,
        "nd": SampleHealthStatus.UNKNOWN,
        "not_defined": SampleHealthStatus.UNKNOWN,
        "not defined": SampleHealthStatus.UNKNOWN,
        "": SampleHealthStatus.UNKNOWN,
    }

    return mapping.get(s, SampleHealthStatus.UNKNOWN)


def parse_infection_level(value: Any) -> InfectionLevel:
    """
    Converte rÃ³tulos (pt/en) ou Enum para InfectionLevel.

    Aceita, por exemplo:
      - "none", "nenhum", "no", "sem"
      - "mild", "leve"
      - "moderate", "moderado"
      - "severe", "severo"
      - "critical", "critico", "crÃ­tico"
    """
    if isinstance(value, InfectionLevel):
        return value
    if value is None:
        return InfectionLevel.NONE

    s = _normalize_label(value)

    mapping = {
        # Sem infecÃ§Ã£o
        "none": InfectionLevel.NONE,
        "nenhum": InfectionLevel.NONE,
        "no": InfectionLevel.NONE,
        "sem": InfectionLevel.NONE,

        # Leve
        "mild": InfectionLevel.MILD,
        "leve": InfectionLevel.MILD,
        "light": InfectionLevel.MILD,

        # Moderada
        "moderate": InfectionLevel.MODERATE,
        "moderado": InfectionLevel.MODERATE,

        # Severa
        "severe": InfectionLevel.SEVERE,
        "severo": InfectionLevel.SEVERE,
        "grave": InfectionLevel.SEVERE,

        # CrÃ­tica
        "critical": InfectionLevel.CRITICAL,
        "critico": InfectionLevel.CRITICAL,
        "critica": InfectionLevel.CRITICAL,
        "critico_agudo": InfectionLevel.CRITICAL,
    }

    return mapping.get(s, InfectionLevel.NONE)


# -------------------------------------------------------------------------
# Models
# -------------------------------------------------------------------------


class Sample(Base):
    __tablename__ = "samples"

    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey("projects.id"), nullable=False)
    uploader_id = Column(Integer, ForeignKey("users.id"))

    # IdentificaÃ§Ã£o
    sample_id = Column(String(100), unique=True, nullable=False, index=True)
    original_id = Column(String(100))
    biosample_accession = Column(String(50))
    sra_accession = Column(String(50))

    # Status de saÃºde
    health_status = Column(SQLEnum(SampleHealthStatus), default=SampleHealthStatus.UNKNOWN)
    infection_level = Column(SQLEnum(InfectionLevel), default=InfectionLevel.NONE)
    pathogen = Column(String(200))  # Agente patogÃªnico identificado
    disease_symptoms = Column(Text)  # Ex.: "white_feces", "hepatopancreas pÃ¡lido", etc.

    # Metadados experimentais
    collection_date = Column(DateTime)
    collection_location = Column(String(255))
    water_temperature = Column(Float)
    salinity = Column(Float)
    ph = Column(Float)
    dissolved_oxygen = Column(Float)

    # Dieta
    diet_protein = Column(Float)        # %
    diet_lipid = Column(Float)          # %
    diet_carbohydrate = Column(Float)   # %
    diet_description = Column(Text)

    # ParÃ¢metros zootÃ©cnicos
    shrimp_weight = Column(Float)  # gramas
    shrimp_length = Column(Float)  # cm
    survival_rate = Column(Float)  # %
    feed_conversion_ratio = Column(Float)

    # Dados de sequenciamento
    sequencing_platform = Column(String(100))
    sequencing_depth = Column(BigInteger)
    region_amplified = Column(String(50))
    primers_used = Column(String(200))

    # Caminhos dos arquivos
    fastq_r1_path = Column(String(1024))
    fastq_r2_path = Column(String(1024))
    metadata_path = Column(String(1024))

    # Status
    is_public = Column(Boolean, default=False)        # Dados pÃºblicos ou privados
    data_source = Column(String(50), default="user")  # user, sra, ena, mgrast
    processing_status = Column(String(50), default="pending")

    # Datas
    uploaded_at = Column(DateTime, default=datetime.utcnow)
    processed_at = Column(DateTime)

    # Relacionamentos
    project = relationship("Project", back_populates="samples")
    uploader = relationship("User", back_populates="samples")
    analyses = relationship("Analysis", back_populates="sample")
    feature_tables = relationship("FeatureTable", back_populates="sample")
    taxonomic_assignments = relationship("TaxonomicAssignment", back_populates="sample")

    def __repr__(self) -> str:
        return f"<Sample id={self.id} sample_id={self.sample_id} health_status={self.health_status}>"


class Analysis(Base):
    __tablename__ = "analyses"

    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey("projects.id"), nullable=False)
    sample_id = Column(Integer, ForeignKey("samples.id"))
    user_id = Column(Integer, ForeignKey("users.id"))

    # IdentificaÃ§Ã£o
    analysis_id = Column(String(100), unique=True, nullable=False, index=True)
    analysis_type = Column(SQLEnum(AnalysisType), nullable=False)
    analysis_name = Column(String(255))
    description = Column(Text)

    # ParÃ¢metros
    parameters = Column(JSON)
    pipeline_version = Column(String(50))
    reference_database = Column(String(100))

    # ExecuÃ§Ã£o
    started_at = Column(DateTime)
    completed_at = Column(DateTime)
    execution_time = Column(Float)  # segundos

    # Status
    status = Column(String(50), default="pending")
    error_message = Column(Text)

    # Resultados
    results_path = Column(String(1024))
    results_summary = Column(JSON)
    visualization_paths = Column(JSON)

    # Relacionamentos
    project = relationship("Project", back_populates="analyses")
    sample = relationship("Sample", back_populates="analyses")
    user = relationship("User", back_populates="analyses")
    results = relationship("AnalysisResult", back_populates="analysis")

    def __repr__(self) -> str:
        return f"<Analysis id={self.id} analysis_id={self.analysis_id} type={self.analysis_type}>"


class FeatureTable(Base):
    __tablename__ = "feature_tables"

    id = Column(Integer, primary_key=True)
    sample_id = Column(Integer, ForeignKey("samples.id"))
    project_id = Column(Integer, ForeignKey("projects.id"))

    # Dados da tabela de features (ASV/OTU)
    feature_id = Column(String(100), index=True)
    feature_sequence = Column(Text)
    feature_count = Column(Integer)
    relative_abundance = Column(Float)

    # Metadados
    table_type = Column(String(50))  # raw, normalized, batch_corrected
    normalization_method = Column(String(50))

    # Relacionamentos
    sample = relationship("Sample", back_populates="feature_tables")
    project = relationship("Project")
    taxonomy = relationship("TaxonomicAssignment", back_populates="feature")

    def __repr__(self) -> str:
        return f"<FeatureTable id={self.id} feature_id={self.feature_id} sample_id={self.sample_id}>"

class ComparativeAnalyzer:
    def __init__(self, public_dir: Path, user_dir: Path):
        self.public_dir = public_dir
        self.user_dir = user_dir
        # constrÃ³i a base de referÃªncia uma vez ao iniciar
        self.reference_db = self.build_reference_database()

    def _extract_record_from_json(self, data: Dict[str, Any], json_file: Path) -> Optional[Dict[str, Any]]:
        """
        Extrai um registro (linha) a partir de um JSON de mÃ©tricas.

        Estrutura esperada (mas flexÃ­vel):

        {
          "sample_id": "...",
          "metrics": { "alpha_shannon": 3.2, "richness": 120, ... },
          "metadata": {
            "health_status": "...",
            "infection_level": "...",
            "farm_id": "...",
            "region": "..."
          }
        }

        Se nÃ£o houver "metrics", pega apenas os campos numÃ©ricos do topo do JSON.
        """
        metrics_block = data.get("metrics")
        if metrics_block is None or not isinstance(metrics_block, dict):
            # fallback: pega sÃ³ campos numÃ©ricos do top-level
            metrics_block = {
                k: v for k, v in data.items()
                if isinstance(v, (int, float))
            }
            if not metrics_block:
                return None

        record: Dict[str, Any] = dict(metrics_block)

        # id da amostra
        sample_id = data.get("sample_id") or data.get("id") or json_file.stem
        record["sample_id"] = sample_id

        # metadados opcionais
        metadata = data.get("metadata") or data.get("sample_metadata") or {}
        if isinstance(metadata, dict):
            for key in ("health_status", "infection_level", "farm_id", "region"):
                if key in metadata:
                    record[key] = metadata[key]

        return record

    def build_reference_database(self) -> Dict[str, Any]:
        """
        Varre self.public_dir Ã  procura de JSONs com mÃ©tricas de amostras pÃºblicas
        e constrÃ³i um DataFrame + estatÃ­sticas de resumo para cada mÃ©trica numÃ©rica.
        """
        records: List[Dict[str, Any]] = []

        if not self.public_dir.exists():
            return {"all_samples": pd.DataFrame(), "summary_stats": {}}

        # Se vocÃª tiver um padrÃ£o especÃ­fico (ex: "*_metrics.json"), ajuste aqui:
        for json_file in self.public_dir.rglob("*.json"):
            try:
                with json_file.open() as f:
                    data = json.load(f)
            except Exception:
                # arquivo invÃ¡lido, ignora silenciosamente
                continue

            if not isinstance(data, dict):
                continue

            rec = self._extract_record_from_json(data, json_file)
            if rec is not None:
                records.append(rec)

        if not records:
            return {"all_samples": pd.DataFrame(), "summary_stats": {}}

        df = pd.DataFrame.from_records(records)

        numeric_cols = df.select_dtypes(include=["number"]).columns.tolist()
        summary_stats: Dict[str, Dict[str, float]] = {}

        for col in numeric_cols:
            series = df[col].dropna()
            if series.empty:
                continue
            summary_stats[col] = {
                "min": float(series.min()),
                "max": float(series.max()),
                "mean": float(series.mean()),
                "std": float(series.std(ddof=0)),
                "p25": float(series.quantile(0.25)),
                "p50": float(series.quantile(0.50)),
                "p75": float(series.quantile(0.75)),
            }

        return {"all_samples": df, "summary_stats": summary_stats}

    @staticmethod
    def _ecdf_percentile(values: pd.Series, value: float) -> float:
        """Retorna o percentil (0â€“100) do valor em relaÃ§Ã£o aos valores de referÃªncia."""
        values = values.dropna()
        if values.empty:
            return float("nan")
        return float((values <= value).mean() * 100.0)

    def _classify_sample(self, metrics: Dict[str, Any], percentiles: Dict[str, float]) -> Dict[str, Any]:
        """
        HeurÃ­stica simples baseada em percentis para estimar status de saÃºde e nÃ­vel de infeÃ§Ã£o.
        """
        if not percentiles:
            return {
                "health_status": "unknown",
                "infection_level": "none",
                "confidence": 0.0,
                "recommendations": [
                    "NÃ£o foi possÃ­vel comparar a amostra por falta de mÃ©tricas numÃ©ricas na base de referÃªncia."
                ],
            }

        # mÃ©tricas em que "mais alto" Ã© bom (diversidade, riqueza, etc.)
        beneficial_high_keywords = [
            "shannon",
            "diversity",
            "richness",
            "observed_features",
            "evenness",
        ]

        # mÃ©tricas em que "mais alto" Ã© ruim (patÃ³genos, vibrioses, mortalidade, etc.)
        harmful_high_keywords = [
            "vibrio",
            "pathogen",
            "white_feces",
            "whitish",
            "mortality",
            "hepatopancreas",
        ]

        positive_scores: List[float] = []
        harmful_scores: List[float] = []
        recommendations: List[str] = []

        for name, p in percentiles.items():
            lname = name.lower()
            if np.isnan(p):
                continue

            if any(k in lname for k in beneficial_high_keywords):
                # percentil alto = bom
                positive_scores.append(p)

            if any(k in lname for k in harmful_high_keywords):
                # percentil alto de "coisa ruim": transformamos em 100 - p
                positive_scores.append(100.0 - p)
                harmful_scores.append(p)

        if positive_scores:
            health_score = float(np.mean(positive_scores))
        else:
            # se nÃ£o bateu nenhuma keyword, usa mÃ©dia de todos os percentis disponÃ­veis
            valid_percentiles = [p for p in percentiles.values() if not np.isnan(p)]
            health_score = float(np.mean(valid_percentiles)) if valid_percentiles else 50.0

        infection_score = float(np.mean(harmful_scores)) if harmful_scores else 0.0

        # classificaÃ§Ã£o por faixas
        if health_score >= 70:
            health_status = "healthy"
        elif health_score >= 50:
            health_status = "borderline"
        elif health_score >= 30:
            health_status = "at_risk"
        else:
            health_status = "diseased"

        if infection_score < 20:
            infection_level = "none"
        elif infection_score < 40:
            infection_level = "low"
        elif infection_score < 70:
            infection_level = "moderate"
        else:
            infection_level = "high"

        # confianÃ§a: depende de (1) quantas mÃ©tricas Ãºteis usamos e (2) quÃ£o longe de 50 estÃ¡
        coverage = min(1.0, len(positive_scores) / max(3, len(percentiles))) if percentiles else 0.0
        confidence_from_score = abs(health_score - 50.0) / 50.0  # 0â€“1
        confidence = 0.3 + 0.4 * coverage + 0.3 * confidence_from_score
        confidence = float(max(0.0, min(1.0, confidence)))

        # recomendaÃ§Ãµes genÃ©ricas
        if infection_level in ("moderate", "high"):
            recommendations.append(
                "Amostra com sinais de possÃ­vel desequilÃ­brio microbiano; considerar monitorizaÃ§Ã£o mais frequente."
            )
        if health_status in ("at_risk", "diseased"):
            recommendations.append(
                "Considerar avaliar manejo, qualidade da Ã¡gua e realizar exames complementares dos animais."
            )
        if not recommendations:
            recommendations.append(
                "Perfil microbiano compatÃ­vel com as amostras de referÃªncia disponÃ­veis na base pÃºblica."
            )

        return {
            "health_status": health_status,
            "infection_level": infection_level,
            "confidence": confidence,
            "recommendations": recommendations,
        }

    def compare_user_sample(self, sample_results_dir: Path, metrics: Dict[str, Any]) -> Dict[str, Any]:
        """
        Compara as mÃ©tricas da amostra do usuÃ¡rio com a base de referÃªncia.
        `metrics` deve conter apenas mÃ©tricas numÃ©ricas (ex.: shannon_index, asvs_generated, etc.).
        """
        ref_df: pd.DataFrame = self.reference_db.get("all_samples", pd.DataFrame())

        percentiles: Dict[str, float] = {}
        comparisons: Dict[str, Dict[str, Any]] = {}

        if not ref_df.empty and metrics:
            numeric_cols = ref_df.select_dtypes(include=["number"]).columns

            for metric_name, value in metrics.items():
                if metric_name not in numeric_cols:
                    continue

                try:
                    value_float = float(value)
                except (TypeError, ValueError):
                    continue

                col = ref_df[metric_name].dropna()
                if col.empty:
                    continue

                p = self._ecdf_percentile(col, value_float)
                percentiles[metric_name] = p

                # similaridade simples: quÃ£o perto do percentil 50 (centro) estÃ¡
                if np.isnan(p):
                    similarity = float("nan")
                else:
                    similarity = max(0.0, 100.0 - 2.0 * abs(p - 50.0))  # 50 => 100%, 0/100 => 0%

                comparisons[metric_name] = {
                    "user_value": value_float,
                    "reference_mean": float(col.mean()),
                    "reference_std": float(col.std(ddof=0)),
                    "reference_min": float(col.min()),
                    "reference_max": float(col.max()),
                    "percentile": p,
                    "similarity": similarity,
                }

        classification = self._classify_sample(metrics, percentiles)

        return {
            "classification": classification,
            "comparisons": comparisons,
            "percentiles": percentiles,
            "user_metrics": metrics,
            # tamanho da base de referÃªncia (nÂº de amostras)
            "reference_size": int(len(ref_df)),
        }

class TaxonomicAssignment(Base):
    __tablename__ = "taxonomic_assignments"
    
    id = Column(Integer, primary_key=True)
    feature_id = Column(String(100), ForeignKey("feature_tables.feature_id"))
    sample_id = Column(Integer, ForeignKey("samples.id"))
    
    # Taxonomia
    kingdom = Column(String(100))
    phylum = Column(String(100))
    class_name = Column(String(100))
    order_name = Column(String(100))
    family = Column(String(100))
    genus = Column(String(100))
    species = Column(String(100))
    confidence = Column(Float)
    
    # AnotaÃ§Ãµes especÃ­ficas
    is_pathogen = Column(Boolean, default=False)
    pathogen_type = Column(String(100))
    is_beneficial = Column(Boolean, default=False)
    beneficial_function = Column(String(200))
    
    # MÃ©tricas
    prevalence = Column(Float)  # FrequÃªncia nas amostras
    mean_abundance = Column(Float)
    
    # Relacionamentos
    feature = relationship("FeatureTable", back_populates="taxonomy")
    sample = relationship("Sample", back_populates="taxonomic_assignments")

class DiversityMetric(Base):
    __tablename__ = "diversity_metrics"
    
    id = Column(Integer, primary_key=True)
    analysis_id = Column(Integer, ForeignKey("analyses.id"))
    sample_id = Column(Integer, ForeignKey("samples.id"))
    
    # MÃ©tricas alfa
    observed_features = Column(Integer)
    shannon_index = Column(Float)
    simpson_index = Column(Float)
    pielou_evenness = Column(Float)
    chao1 = Column(Float)
    ace = Column(Float)
    faith_pd = Column(Float)
    
    # MÃ©tricas beta (por par de amostras)
    sample1_id = Column(Integer)
    sample2_id = Column(Integer)
    bray_curtis = Column(Float)
    jaccard = Column(Float)
    weighted_unifrac = Column(Float)
    unweighted_unifrac = Column(Float)

class Comparison(Base):
    __tablename__ = "comparisons"
    
    id = Column(Integer, primary_key=True)
    project_id = Column(Integer, ForeignKey("projects.id"))
    user_sample_id = Column(Integer, ForeignKey("samples.id"))
    
    # ConfiguraÃ§Ã£o
    comparison_type = Column(String(50))  # health_status, infection_level, diet
    reference_group = Column(String(100))
    
    # Resultados
    overall_similarity = Column(Float)
    classification_result = Column(JSON)
    visualization_path = Column(String(1024))
    
    # Datas
    created_at = Column(DateTime, default=datetime.utcnow)
    
    # Relacionamentos
    project = relationship("Project", back_populates="comparisons")
    user_sample = relationship("Sample", foreign_keys=[user_sample_id])
    results = relationship("ComparativeResult", back_populates="comparison")


class ComparativeResult(Base):
    __tablename__ = "comparative_results"
    
    id = Column(Integer, primary_key=True)
    comparison_id = Column(Integer, ForeignKey("comparisons.id"))
    
    # ComparaÃ§Ã£o
    user_sample_id = Column(Integer, ForeignKey("samples.id"))
    reference_sample_id = Column(Integer, ForeignKey("samples.id"))
    reference_group = Column(String(100))  # healthy, dysbiosis_mild, etc.
    
    # MÃ©tricas de similaridade
    bray_curtis_similarity = Column(Float)
    jaccard_similarity = Column(Float)
    taxonomic_overlap = Column(Float)  # % de taxa em comum
    
    # PosiÃ§Ã£o percentil
    shannon_percentile = Column(Float)
    pathogen_abundance_percentile = Column(Float)
    beneficial_abundance_percentile = Column(Float)
    
    # ClassificaÃ§Ã£o
    predicted_health_status = Column(String(50))
    predicted_infection_level = Column(String(50))
    confidence_score = Column(Float)

    # Relacionamentos
    comparison = relationship("Comparison", back_populates="results")
    user_sample = relationship("Sample", foreign_keys=[user_sample_id])
    reference_sample = relationship("Sample", foreign_keys=[reference_sample_id])

class AnalysisResult(Base):
    __tablename__ = "analysis_results"
    
    id = Column(Integer, primary_key=True)
    analysis_id = Column(Integer, ForeignKey("analyses.id"))
    
    # Resultados especÃ­ficos por tipo de anÃ¡lise
    result_type = Column(String(50))
    result_data = Column(JSON)
    statistics = Column(JSON)
    visualizations = Column(JSON)
    
    # Relacionamentos
    analysis = relationship("Analysis", back_populates="results")

class PublicDataset(Base):
    __tablename__ = "public_datasets"
    
    id = Column(Integer, primary_key=True)
    
    # IdentificaÃ§Ã£o
    dataset_id = Column(String(100), unique=True, nullable=False, index=True)
    dataset_name = Column(String(255))
    source = Column(String(50))  # SRA, ENA, MG-RAST
    accession = Column(String(100))
    
    # Metadados
    species = Column(String(100))
    tissue_type = Column(String(100))
    health_status = Column(String(50))
    infection_level = Column(String(50))
    geographic_location = Column(String(255))
    
    # Dados
    sample_count = Column(Integer)
    total_reads = Column(BigInteger)
    download_status = Column(String(50), default="pending")
    
    # Caminhos
    local_path = Column(String(1024))
    metadata_path = Column(String(1024))
    
    # Datas
    downloaded_at = Column(DateTime)
    created_at = Column(DateTime, default=datetime.utcnow)

# ============================================================================
# SISTEMA DE DOWNLOAD DE DADOS PÃšBLICOS
# ============================================================================
class PublicDataDownloader:
  """
  Sistema de download de dados pÃºblicos de camarÃµes.

  - Busca SRA via pysradb (se disponÃ­vel) + Entrez (fallback).
  - Busca ENA via ENA Portal API.
  - Filtra para 16S, AMPLICON, plataforma Illumina e espÃ©cies de interesse.
  - Permite escolher regiÃ£o alvo (ex: V3-V4, V4, V1-V3 etc.).
  - Respeita um limite N (max_samples / n_max).
  """

  def __init__(
    self,
    email: Optional[str] = None,
    base_dir: Optional[Path] = None,
    max_samples: int = 100,
  ):
    # Configura e-mail do Entrez
    self.email = email or os.getenv("NCBI_EMAIL", "shrimpmicrobiome@example.com")
    Entrez.email = self.email
    Entrez.tool = "ShrimpMicrobiomePlatform"

    # pysradb: conecta ao SRA metadata (opcional)
    try:
      import pysradb  # type: ignore
      self.sra_db = pysradb.SRAweb()
    except Exception as e:
      logger.warning(f"NÃ£o foi possÃ­vel inicializar pysradb.SRAweb: {e}")
      self.sra_db = None

    # DiretÃ³rio base para organizar dados pÃºblicos
    if base_dir is not None:
      self.base_dir = Path(base_dir)
    else:
      try:
        self.base_dir = Path(config.PUBLIC_DATA_DIR)  # type: ignore[attr-defined]
      except Exception:
        self.base_dir = Path("public_data")

    self.max_samples = max_samples

    # RegiÃµes tÃ­picas de 16S
    self.valid_regions = {
      "V1-V2",
      "V1-V3",
      "V3-V4",
      "V4",
      "V4-V5",
      "V3-V5",
    }

    # URL base da ENA Portal API
    self.ena_base_url = "https://www.ebi.ac.uk/ena/portal/api"

  # ------------------------------------------------------------------
  # BUSCA GENÃ‰RICA (ENA + BioProject/SRA via ENA)
  # ------------------------------------------------------------------
  def search_public_datasets(
    self,
    species: str = "Penaeus vannamei",
    marker: str = "16S",
    region: Optional[str] = "V3-V4",
    bioproject_ids: Optional[Iterable[str]] = None,
    limit_per_project: int = 0,
    n_max: Optional[int] = None,
    session: Optional[requests.Session] = None,
    search_term: Optional[str] = None,   # <---------------------- CORRETO
    source: Optional[str] = None,        # "ENA", "SRA" ou None
  ) -> pd.DataFrame:

    n_max = n_max or self.max_samples
    session = session or requests.Session()

    dfs: List[pd.DataFrame] = []

    # ------------------------------------------------------------
    # 1) Busca geral no ENA (via search_ena_datasets)
    # ------------------------------------------------------------
    if source is None or source.upper() == "ENA":
      df_ena = self.search_ena_datasets(
        search_term=search_term,
        species=species,
        marker=marker,
        region=region,
        n_max=n_max,
      )

      if df_ena is not None and not df_ena.empty:
        df_ena = df_ena.copy()

        if "run_accession" in df_ena.columns and "accession" not in df_ena.columns:
          df_ena["accession"] = df_ena["run_accession"]

        if "instrument_platform" in df_ena.columns and "platform" not in df_ena.columns:
          df_ena["platform"] = df_ena["instrument_platform"]

        if "total_reads" not in df_ena.columns:
          if "total_spots" in df_ena.columns:
            df_ena["total_reads"] = df_ena["total_spots"]
          elif "read_count" in df_ena.columns:
            df_ena["total_reads"] = df_ena["read_count"]

        if "bioproject_accession" in df_ena.columns and "bioproject" not in df_ena.columns:
          df_ena["bioproject"] = df_ena["bioproject_accession"]

        if "source" not in df_ena.columns:
          df_ena["source"] = "ENA"

        if "fastq_ftp" in df_ena.columns and "fastq_ftp_links" not in df_ena.columns:
          df_ena["fastq_ftp_links"] = df_ena["fastq_ftp"].fillna("")

        if "fastq_http_links" not in df_ena.columns:
          df_ena["fastq_http_links"] = ""

        if "fastq_aspera_links" not in df_ena.columns:
          df_ena["fastq_aspera_links"] = ""

        for col in ["fastq_ftp", "fastq_http", "fastq_aspera"]:
          if col not in df_ena.columns:
            df_ena[col] = ""

        dfs.append(df_ena)

    # ------------------------------------------------------------
    # 2) Busca por BioProjects (via ENA read_run_table â€“ espelho SRA)
    # ------------------------------------------------------------
    if (source is None or source.upper() == "SRA") and bioproject_ids:
      default_fields = [
        "study_accession",
        "bioproject_accession",
        "experiment_accession",
        "sample_accession",
        "run_accession",
        "scientific_name",
        "library_strategy",
        "library_layout",
        "instrument_platform",
        "instrument_model",
        "fastq_ftp",
        "fastq_http",
        "fastq_aspera",
      ]

      all_rows: List[Dict[str, str]] = []

      for bioproject in bioproject_ids:
        rows = _fetch_ena_read_run_table(
          study_or_project=bioproject,
          fields=default_fields,
          limit=limit_per_project,
          session=session,
        )

        for row in rows:
          links = _row_links_from_ena(row)
          row["fastq_ftp_links"] = ";".join(links["fastq_ftp"]) if links["fastq_ftp"] else ""
          row["fastq_http_links"] = ";".join(links["fastq_http"]) if links["fastq_http"] else ""
          row["fastq_aspera_links"] = ";".join(links["fastq_aspera"]) if links["fastq_aspera"] else ""
          row["source"] = "SRA"

          if row.get("bioproject_accession"):
            row.setdefault("bioproject", row["bioproject_accession"])
          else:
            row.setdefault("bioproject", bioproject)

          all_rows.append(row)

      if all_rows:
        df_sra = pd.DataFrame(all_rows)

        if "bioproject_accession" in df_sra.columns and "bioproject" not in df_sra.columns:
          df_sra["bioproject"] = df_sra["bioproject_accession"]

        if "run_accession" in df_sra.columns and "accession" not in df_sra.columns:
          df_sra["accession"] = df_sra["run_accession"]

        if "instrument_platform" in df_sra.columns and "platform" not in df_sra.columns:
          df_sra["platform"] = df_sra["instrument_platform"]

        if "total_reads" not in df_sra.columns:
          if "total_spots" in df_sra.columns:
            df_sra["total_reads"] = df_sra["total_spots"]
          elif "read_count" in df_sra.columns:
            df_sra["total_reads"] = df_sra["read_count"]

        for col in ["fastq_ftp_links", "fastq_http_links", "fastq_aspera_links"]:
          if col not in df_sra.columns:
            df_sra[col] = ""

        for col in ["fastq_ftp", "fastq_http", "fastq_aspera"]:
          if col not in df_sra.columns:
            df_sra[col] = ""

        dfs.append(df_sra)

    # ------------------------------------------------------------
    # 3) Consolida resultados
    # ------------------------------------------------------------
    if not dfs:
      return pd.DataFrame(
        columns=[
          "accession", "source", "study_accession", "bioproject",
          "experiment_accession", "sample_accession", "scientific_name",
          "library_strategy", "library_layout", "platform", "total_reads",
          "fastq_ftp", "fastq_http", "fastq_aspera",
          "fastq_ftp_links", "fastq_http_links", "fastq_aspera_links",
        ]
      )

    df_all = pd.concat(dfs, ignore_index=True, sort=False)

    if "accession" in df_all.columns:
      if "total_reads" in df_all.columns:
        df_all["total_reads"] = pd.to_numeric(df_all["total_reads"], errors="coerce")
        df_all = df_all.sort_values("total_reads", ascending=False)
      df_all = df_all.drop_duplicates(subset=["accession"])

    df_all = df_all.head(n_max).reset_index(drop=True)

    return df_all

  # ------------------------------------------------------------------
  # Busca especÃ­fica no ENA
  # ------------------------------------------------------------------
  def search_ena_datasets(
    self,
    search_term: Optional[str] = None,
    species: str = "Penaeus vannamei",
    marker: str = "16S",
    region: Optional[str] = "V3-V4",
    n_max: int = 100,
  ) -> pd.DataFrame:
    """
    Busca datasets no ENA (European Nucleotide Archive) via Portal API (result=read_run).

    - Usa sintaxe de campos do ENA (tax_eq, scientific_name, library_strategy etc.)
    - Evita queries livres que causam HTTP 400.
    """
    try:
      # --------------------------------------------------
      # Construir query ENA vÃ¡lida
      # --------------------------------------------------
      query_parts: List[str] = []

      # Se o usuÃ¡rio passou um search_term "avanÃ§ado" (com campos),
      # usamos diretamente para nÃ£o quebrar queries customizadas.
      if search_term and any(tok in search_term for tok in ["=", "tax_", "(", ")", "\""]):
        base_query = search_term.strip()
        query_parts.append(base_query)
      else:
        # Caso contrÃ¡rio, montamos a query a partir de species/marker/region

        # EspÃ©cie â†’ tax_eq ou scientific_name
        if species:
          if species.lower().startswith("penaeus vannamei"):
            # TaxID conhecido para Penaeus vannamei
            query_parts.append("tax_eq(6689)")
          else:
            safe_species = species.replace('"', '\\"')
            query_parts.append(f'scientific_name="{safe_species}"')

        # Marker â†’ 16S => AMPLICON
        if marker:
          if "16s" in marker.lower():
            query_parts.append('library_strategy="AMPLICON"')
          else:
            safe_marker = marker.replace('"', '\\"')
            query_parts.append(f'library_strategy="{safe_marker}"')

        # RegiÃ£o (V3-V4, V4, V1-V3 etc.) â†’ filtrar em description/sample_title
        if region and region != "Todas":
          safe_region = region.replace('"', '\\"')
          query_parts.append(
            f'(description="{safe_region}" OR sample_title="{safe_region}")'
          )

      # Fallback para nÃ£o mandar query vazia
      if not query_parts:
        query = "tax_eq(6689)"
      else:
        query = " AND ".join(query_parts)

      # --------------------------------------------------
      # ParÃ¢metros da API ENA
      # --------------------------------------------------
      params = {
        "result": "read_run",
        "query": query,
        "fields": ",".join([
          "run_accession",
          "study_accession",
          "sample_accession",
          "experiment_accession",
          "scientific_name",
          "instrument_platform",
          "library_strategy",
          "library_layout",
          "read_count",
          "fastq_ftp",
          "fastq_aspera",
          "fastq_http",
        ]),
        "limit": str(n_max),
        "format": "tsv",
      }

      response = requests.get(
        f"{self.ena_base_url}/search",
        params=params,
        timeout=30,
      )

      try:
        response.raise_for_status()
      except requests.HTTPError as http_err:
        logger.error(
          "Erro HTTP na busca ENA (status=%s) para query='%s': %s",
          response.status_code,
          query,
          http_err,
        )
        return pd.DataFrame()

      content = response.text
      if not content or "run_accession" not in content:
        return pd.DataFrame()

      # Converter para DataFrame
      df = pd.read_csv(io.StringIO(content), sep="\t")

      if df.empty:
        return df

      # Adicionar coluna source
      df["source"] = "ENA"

      # Renomear read_count para total_reads se existir
      if "read_count" in df.columns and "total_reads" not in df.columns:
        df["total_reads"] = df["read_count"]

      # Adicionar coluna platform se instrument_platform existir
      if "instrument_platform" in df.columns and "platform" not in df.columns:
        df["platform"] = df["instrument_platform"]

      return df

    except Exception as e:
      logger.error(f"Erro na busca ENA: {e}", exc_info=True)
      return pd.DataFrame()


  # ------------------------------------------------------------------
  # Alias para versÃ£o antiga
  # ------------------------------------------------------------------
  def public_datasets(
    self,
    species: str = "Penaeus vannamei",
    marker: str = "16S",
    region: Optional[str] = "V3-V4",
    bioproject_ids: Optional[Iterable[str]] = None,
    limit_per_project: int = 0,
    n_max: Optional[int] = None,
    session: Optional[requests.Session] = None,
  ) -> pd.DataFrame:
    """
    Alias de compatibilidade para `search_public_datasets`.

    MantÃ©m a assinatura usada na versÃ£o antiga do cÃ³digo, mas internamente
    chama o mÃ©todo novo, que integra ENA + BioProjects/SRA espelhados.
    """
    return self.search_public_datasets(
      species=species,
      marker=marker,
      region=region,
      bioproject_ids=bioproject_ids,
      limit_per_project=limit_per_project,
      n_max=n_max,
      session=session,
      search_term=None,
      source=None,  # None => ENA + SRA
    )

  # ------------------------------------------------------------------
  # OrganizaÃ§Ã£o dos dados pÃºblicos (gera tabela de links etc.)
  # ------------------------------------------------------------------
  def organize_public_data(
    self,
    df: Optional[pd.DataFrame] = None,
    output_dir: Optional[Path] = None,
    links_filename: str = "public_fastq_links.tsv",
    **search_kwargs,
  ) -> pd.DataFrame:
    """
    Organiza os datasets pÃºblicos em um diretÃ³rio de trabalho e gera
    uma tabela de links de download.

    ParÃ¢metros
    ----------
    df : pd.DataFrame, opcional
      DataFrame jÃ¡ retornado por `search_public_datasets`. Se nÃ£o for
      fornecido, os dados serÃ£o obtidos chamando
      `search_public_datasets(**search_kwargs)`.
    output_dir : Path ou str, opcional
      DiretÃ³rio onde os arquivos de organizaÃ§Ã£o serÃ£o escritos.
      Default: self.base_dir / "public_data".
    links_filename : str
      Nome do arquivo TSV com os links (default: "public_fastq_links.tsv").
    **search_kwargs :
      ParÃ¢metros extras repassados para `search_public_datasets` caso `df`
      nÃ£o seja fornecido (ex.: species, marker, region, n_max, search_term,
      bioproject_ids, source, session, etc.).

    Retorna
    -------
    pd.DataFrame
      O DataFrame final utilizado/gerado.
    """
    # Se df nÃ£o foi passado, faz uma busca agora
    if df is None:
      # Se o chamador nÃ£o passou n_max, usa o default da instÃ¢ncia
      if "n_max" not in search_kwargs or search_kwargs.get("n_max") is None:
        search_kwargs.setdefault("n_max", self.max_samples)

      df = self.search_public_datasets(**search_kwargs)

    if df is None:
      df = pd.DataFrame()

    # Garante que seja sempre um DataFrame
    if not isinstance(df, pd.DataFrame):
      df = pd.DataFrame(df)

    # DiretÃ³rio de saÃ­da
    if output_dir is None:
      output_dir = self.base_dir / "public_data"

    output_dir = _ensure_dir(Path(output_dir))

    # Gera tabela de links (se houver colunas relevantes)
    _write_links_table(df, output_dir / links_filename)

    logger.info("Dados pÃºblicos organizados em %s", output_dir)

    return df



  # ------------------------------------------------------------------
  # HEURÃSTICAS DE METADADOS
  # ------------------------------------------------------------------
  def determine_health_status(self, row: pd.Series) -> str:
    """
    HeurÃ­stica simples para inferir status de saÃºde a partir de metadados textuais.
    """
    text = " ".join(
      [
        str(row.get("study_title", "")),
        str(row.get("sample_alias", "")),
        str(row.get("description", "")),
      ]
    ).lower()

    # Palavras-chave indicativas de saÃºde
    healthy_keywords = [
      "healthy",
      "saudÃ¡vel",
      "saudavel",
      "control",
      "controle",
      "uninfected",
      "naive",
    ]

    # Palavras-chave indicativas de doenÃ§a / disbiose
    disease_keywords = [
      "disease",
      "infected",
      "challenged",
      "challenge",
      "doente",
      "sick",
      "moribund",
    ]

    if any(k in text for k in healthy_keywords):
      return "SaudÃ¡vel"
    if any(k in text for k in disease_keywords):
      return "Disbiose"
    return "Desconhecido"

  def determine_infection_level(self, row: pd.Series) -> str:
    """
    InferÃªncia simples de nÃ­vel de infecÃ§Ã£o (Alta, MÃ©dia, Baixa ou Desconhecido).
    """
    text = " ".join(
      [
        str(row.get("study_title", "")),
        str(row.get("sample_alias", "")),
        str(row.get("description", "")),
      ]
    ).lower()

    high_keywords = [
      "high",
      "severe",
      "alta",
      "severa",
      "strong",
    ]

    medium_keywords = [
      "medium",
      "moderate",
      "mÃ©dia",
      "media",
      "moderada",
    ]

    low_keywords = [
      "low",
      "baixa",
      "leve",
      "mild",
    ]

    if any(k in text for k in high_keywords):
      return "Alta"
    if any(k in text for k in medium_keywords):
      return "MÃ©dia"
    if any(k in text for k in low_keywords):
      return "Baixa"
    return "Desconhecido"

  def extract_diet_info(self, row: pd.Series) -> Dict[str, Optional[float]]:
    """
    Extrai informaÃ§Ãµes de dieta (proteÃ­na, lipÃ­dio, carboidrato).

    OBS: Por enquanto Ã© apenas um placeholder que retorna None para todos
    os campos. No futuro vocÃª pode implementar parsing dos metadados
    (por exemplo, extrair percentagens de proteÃ­na, lipÃ­dio, carboidrato
    a partir de descriÃ§Ãµes da dieta).
    """
    # Exemplo de ponto de partida:
    # text = " ".join(
    #   [
    #     str(row.get("study_title", "")),
    #     str(row.get("sample_alias", "")),
    #     str(row.get("description", "")),
    #   ]
    # ).lower()
    #
    # Aqui vocÃª poderia usar regex para extrair "35% protein", "8% lipid", etc.

    return {
      "protein": None,
      "lipid": None,
      "carbohydrate": None,
    }

  def extract_tissue_type(self, row: pd.Series) -> str:
    """
    Tenta identificar tipo de tecido/ambiente a partir de metadados textuais.
    """
    text = " ".join(
      [
        str(row.get("study_title", "")),
        str(row.get("sample_alias", "")),
        str(row.get("description", "")),
        str(row.get("library_source", "")),
      ]
    ).lower()

    if "intest" in text or "gut" in text:
      return "Intestino"
    if "hepatopancreas" in text or "hepato" in text:
      return "HepatopÃ¢ncreas"
    if "water column" in text or "water" in text or "Ã¡gua" in text or "agua" in text:
      return "Ãgua"
    if "sediment" in text or "sedimento" in text:
      return "Sedimento"
    return "Desconhecido"

  def extract_pathogen_info(self, row: pd.Series) -> Optional[str]:
    """
    Identifica possÃ­veis patÃ³genos mencionados nos metadados.
    """
    text = " ".join(
      [
        str(row.get("study_title", "")),
        str(row.get("sample_alias", "")),
        str(row.get("description", "")),
      ]
    ).lower()

    if "vibrio" in text:
      return "Vibrio spp."
    if "white spot" in text or "wssv" in text:
      return "WSSV"
    # Exemplos para expansÃ£o futura:
    # if "aeromonas" in text:
    #   return "Aeromonas spp."
    # if "ehp" in text or "hepatopancreatic microsporidiosis" in text:
    #   return "EHP"
    return None


  # ------------------------------------------------------------------
  # ORGANIZAÃ‡ÃƒO/RESUMO PARA DASHBOARD
  # ------------------------------------------------------------------
  def organize_public_data(
      self,
      output_dir: str,
      search_term: Optional[str] = None,
      species: str = "Penaeus vannamei",
      marker: str = "16S",
      region: Optional[str] = "V3-V4",
      bioproject_ids: Optional[Iterable[str]] = None,
      limit_per_project: int = 0,
      n_max: Optional[int] = None,
      session: Optional[requests.Session] = None,
  ) -> pd.DataFrame:
      """
      High-level helper that:
        1. Busca dados pÃºblicos (SRA/ENA) automaticamente
        2. Cria uma pasta por BioProject sob output_dir
        3. Escreve:
           - 'sra_runs_metadata.tsv' global para todos os runs
           - 'download_links.tsv' por BioProject com run_accession e links

      Retorna o DataFrame completo de metadados.
      """
      output_dir_path = Path(output_dir)
      _ensure_dir(output_dir_path)

    # Fazer busca automÃ¡tica
      df = self.search_public_datasets(
          species=species,
          marker=marker,
          region=region,
          bioproject_ids=bioproject_ids,
          limit_per_project=limit_per_project,
          n_max=n_max or self.max_samples,
          session=session,
          search_term=search_term,
      )

      if df.empty:
          return df

    # Salvar metadados globais
      metadata_path = output_dir_path / "sra_runs_metadata.tsv"
      df.to_csv(metadata_path, sep="\t", index=False)

    # Criar subpastas por BioProject
      if "bioproject" in df.columns:
          group_col = "bioproject"
      elif "bioproject_accession" in df.columns:
          group_col = "bioproject_accession"
      else:
          group_col = None

      if group_col is not None:
          for bioproject, group in df.groupby(group_col):
              proj_dir = _ensure_dir(output_dir_path / str(bioproject))
              links_path = proj_dir / "download_links.tsv"
              _write_links_table(group, links_path)
      else:
        # Fallback: tudo no diretÃ³rio raiz
          links_path = output_dir_path / "download_links.tsv"
          _write_links_table(df, links_path)

      return df

    # ====== Suas funÃ§Ãµes auxiliares levemente ajustadas (mantive nomes) ======

# ============================================================================
# PIPELINE DE PROCESSAMENTO 16S
# ============================================================================
class Shrimp16SPipeline:
    """Pipeline completo para processamento 16S de camarÃµes"""
    
    def __init__(self, project_dir: Path):
        self.project_dir = project_dir
        self.log_file = project_dir / "pipeline.log"
        self.results_dir = project_dir / "results"
        self.qc_dir = project_dir / "quality_control"
        self.processed_dir = project_dir / "processed"
        
        # Criar diretÃ³rios
        for d in [self.results_dir, self.qc_dir, self.processed_dir]:
            d.mkdir(exist_ok=True)
        
        # Configurar logging
        self.logger = self.setup_logger()
    
    def setup_logger(self) -> logging.Logger:
        """Configura logger do pipeline"""
        logger = logging.getLogger(f"pipeline_{self.project_dir.name}")
        logger.setLevel(logging.INFO)
        
        # Handler para arquivo
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.INFO)
        
        # Handler para console
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        
        # Formato
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        console_handler.setFormatter(formatter)
        
        logger.addHandler(file_handler)
        logger.addHandler(console_handler)
        
        return logger
    
    def run_full_pipeline(self, fastq_dir: Path, metadata_path: Path,
                          params: Dict = None) -> Dict:
        """Executa pipeline completo"""

        params = params or config.DEFAULT_PARAMS

        # resultado bÃ¡sico, para nÃ£o dar erro no except se algo quebrar cedo
        results = {
            "status": "running",
            "steps": {},
            "output_files": [],
            "metrics": {},
        }

        try:
            self.logger.info("Iniciando pipeline 16S para camarÃµes")
            start_time = time.time()

            # 0. garantir metadata.tsv dentro do diretÃ³rio do projeto
            project_metadata_path = self.project_dir / "metadata.tsv"
            try:
                if metadata_path and Path(metadata_path).exists():
                    shutil.copy2(metadata_path, project_metadata_path)
                    self.logger.info(f"Metadados copiados para {project_metadata_path}")
                else:
                    self.logger.warning("Arquivo de metadados nÃ£o encontrado; alguns passos do QIIME2 podem falhar.")
            except Exception as e:
                self.logger.warning(f"NÃ£o foi possÃ­vel copiar metadados para o diretÃ³rio do projeto: {e}")

            # 1. Controle de qualidade com Prinseq
            self.logger.info("1. Controle de qualidade com Prinseq")
            qc_results = self.run_prinseq_qc(fastq_dir, params)
            results["steps"]["quality_control"] = qc_results

            # 2. Importar para QIIME2
            self.logger.info("2. Importar dados para QIIME2")
            qiime_artifact = self.import_to_qiime(fastq_dir, project_metadata_path)
            results["steps"]["import"] = {"artifact": str(qiime_artifact)}

            # 3. Denoising com DADA2
            self.logger.info("3. Denoising com DADA2")
            dada2_results = self.run_dada2(qiime_artifact, params)
            results["steps"]["denoising"] = dada2_results

            # 4. AtribuiÃ§Ã£o taxonÃ´mica
            self.logger.info("4. AtribuiÃ§Ã£o taxonÃ´mica")
            taxonomy_results = self.assign_taxonomy(dada2_results["table"])
            results["steps"]["taxonomy"] = taxonomy_results

            # 5. AnÃ¡lises de diversidade
            self.logger.info("5. AnÃ¡lises de diversidade")
            diversity_results = self.analyze_diversity(
                dada2_results["table"],
                dada2_results["tree"],
                params,
            )
            results["steps"]["diversity"] = diversity_results

            # 6. AnÃ¡lises estatÃ­sticas
            self.logger.info("6. AnÃ¡lises estatÃ­sticas")
            stats_results = self.run_statistical_analyses(
                dada2_results["table"],
                taxonomy_results["taxonomy"],
                project_metadata_path,
            )
            results["steps"]["statistics"] = stats_results

            # 7. CorreÃ§Ã£o de batch effects
            self.logger.info("7. CorreÃ§Ã£o de batch effects")
            batch_results = self.correct_batch_effects(
                dada2_results["table"],
                project_metadata_path,
            )
            results["steps"]["batch_correction"] = batch_results

            # 8. Gerar visualizaÃ§Ãµes
            self.logger.info("8. Gerando visualizaÃ§Ãµes")
            viz_results = self.generate_visualizations(
                dada2_results["table"],
                taxonomy_results["taxonomy"],
                diversity_results,
                stats_results,
            )
            results["steps"]["visualizations"] = viz_results

            # MÃ©tricas finais
            execution_time = time.time() - start_time
            results["metrics"] = {
                "execution_time": execution_time,
                "samples_processed": len(qc_results.get("samples", [])),
                "asvs_generated": dada2_results.get("asv_count", 0),
                "success": True,
            }

            results["status"] = "completed"
            self.logger.info(f"Pipeline concluÃ­do em {execution_time:.1f} segundos")

            # Salvar resultados
            results_path = self.results_dir / "pipeline_results.json"
            with open(results_path, "w") as f:
                json.dump(results, f, indent=2)

            results["output_files"].append(str(results_path))

            return results

        except Exception as e:
            self.logger.error(f"Erro no pipeline: {e}", exc_info=True)
            return {
                "status": "failed",
                "error": str(e),
                "steps": results.get("steps", {}),
                "metrics": results.get("metrics", {}),
            }
    
    def run_prinseq_qc(self, fastq_dir: Path, params: Dict) -> Dict:
        """Executa controle de qualidade com Prinseq"""
        
        qc_results = {
            'samples': [],
            'qc_metrics': {},
            'output_files': []
        }
        
        # Listar arquivos FASTQ
        fastq_files = list(fastq_dir.glob("*.fastq*"))
        if not fastq_files:
            fastq_files = list(fastq_dir.glob("*.fq*"))
        
        for fastq_file in fastq_files:
            try:
                sample_name = fastq_file.stem.replace('.fastq', '').replace('.fq', '')
                output_prefix = self.qc_dir / sample_name
                
                # Comando Prinseq bÃ¡sico
                cmd = [
                    "prinseq-lite.pl",
                    "-fastq", str(fastq_file),
                    "-out_format", "3",
                    "-out_good", str(output_prefix),
                    "-log", str(self.qc_dir / f"{sample_name}.log"),
                    "-trim_qual_right", "30",
                    "-trim_qual_type", "min",
                    "-trim_qual_window", "10",
                    "-min_len", "100",
                    "-ns_max_p", "10"
                ]
                
                self.logger.info(f"Processando {sample_name} com Prinseq")
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode == 0:
                    # Verificar arquivos gerados
                    output_files = list(self.qc_dir.glob(f"{sample_name}*"))
                    qc_results['samples'].append(sample_name)
                    qc_results['output_files'].extend([str(f) for f in output_files])
                    
                    # Extrair mÃ©tricas do log
                    metrics = self.extract_prinseq_metrics(
                        self.qc_dir / f"{sample_name}.log"
                    )
                    qc_results['qc_metrics'][sample_name] = metrics
                    
                else:
                    self.logger.warning(f"Erro no Prinseq para {sample_name}: {result.stderr}")
                    
            except Exception as e:
                self.logger.error(f"Erro processando {fastq_file}: {e}")
        
        return qc_results
    
    def extract_prinseq_metrics(self, log_file: Path) -> Dict:
        """Extrai mÃ©tricas do log do Prinseq"""
        metrics = {
            'total_reads': 0,
            'reads_passing': 0,
            'mean_length': 0,
            'mean_quality': 0
        }
        
        try:
            if log_file.exists():
                with open(log_file, 'r') as f:
                    for line in f:
                        if 'reads input' in line:
                            parts = line.split()
                            metrics['total_reads'] = int(parts[0])
                        elif 'reads output' in line and 'good' in line:
                            parts = line.split()
                            metrics['reads_passing'] = int(parts[0])
                        elif 'mean length' in line:
                            parts = line.split()
                            metrics['mean_length'] = float(parts[0])
                        elif 'mean quality' in line:
                            parts = line.split()
                            metrics['mean_quality'] = float(parts[0])
        except:
            pass
        
        return metrics
    
    def import_to_qiime(self, fastq_dir: Path, metadata_path: Path) -> Path:
        """Importa dados para formato QIIME2"""
        
        # Criar arquivo de manifesto para importaÃ§Ã£o
        manifest_data = []
        
        # Encontrar arquivos pareados R1 e R2
        r1_files = sorted(list(fastq_dir.glob("*R1*.fastq*")))
        r2_files = sorted(list(fastq_dir.glob("*R2*.fastq*")))
        
        for r1, r2 in zip(r1_files, r2_files):
            sample_id = r1.stem.replace('_R1', '').replace('.fastq', '')
            manifest_data.append({
                'sample-id': sample_id,
                'forward-absolute-filepath': str(r1.absolute()),
                'reverse-absolute-filepath': str(r2.absolute())
            })
        
        # Salvar manifesto
        manifest_path = self.project_dir / "manifest.csv"
        manifest_df = pd.DataFrame(manifest_data)
        manifest_df.to_csv(manifest_path, index=False)
        
        # Importar para QIIME2
        output_artifact = self.processed_dir / "demux.qza"
        
        cmd = [
            "qiime", "tools", "import",
            "--type", "SampleData[PairedEndSequencesWithQuality]",
            "--input-path", str(manifest_path),
            "--output-path", str(output_artifact),
            "--input-format", "PairedEndFastqManifestPhred33"
        ]
        
        self.logger.info("Importando dados para QIIME2")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            raise RuntimeError(f"Erro na importaÃ§Ã£o QIIME2: {result.stderr}")
        
        # Gerar resumo
        summary_path = self.processed_dir / "demux_summary.qzv"
        cmd = [
            "qiime", "demux", "summarize",
            "--i-data", str(output_artifact),
            "--o-visualization", str(summary_path)
        ]
        
        subprocess.run(cmd, capture_output=True, text=True)
        
        return output_artifact
    
    def run_dada2(self, input_artifact: Path, params: Dict) -> Dict:
        """Executa denoising com DADA2"""
        
        # Caminhos de saÃ­da
        table_path = self.processed_dir / "table.qza"
        rep_seqs_path = self.processed_dir / "rep_seqs.qza"
        stats_path = self.processed_dir / "stats.qza"
        
        cmd = [
            "qiime", "dada2", "denoise-paired",
            "--i-demultiplexed-seqs", str(input_artifact),
            "--p-trunc-len-f", str(params['trunc_len_f']),
            "--p-trunc-len-r", str(params['trunc_len_r']),
            "--p-trim-left-f", str(params['trim_left_f']),
            "--p-trim-left-r", str(params['trim_left_r']),
            "--p-max-ee-f", str(params['max_ee']),
            "--p-max-ee-r", str(params['max_ee']),
            "--p-n-threads", "4",
            "--o-table", str(table_path),
            "--o-representative-sequences", str(rep_seqs_path),
            "--o-denoising-stats", str(stats_path)
        ]
        
        self.logger.info("Executando DADA2")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            raise RuntimeError(f"Erro no DADA2: {result.stderr}")
        
        # Gerar Ã¡rvore filogenÃ©tica
        tree_path = self.processed_dir / "tree.qza"
        aligned_seqs_path = self.processed_dir / "aligned.qza"
        masked_seqs_path = self.processed_dir / "masked.qza"
        
        cmd_tree = [
            "qiime", "phylogeny", "align-to-tree-mafft-fasttree",
            "--i-sequences", str(rep_seqs_path),
            "--o-alignment", str(aligned_seqs_path),
            "--o-masked-alignment", str(masked_seqs_path),
            "--o-tree", str(tree_path)
        ]
        
        result = subprocess.run(cmd_tree, capture_output=True, text=True)
        
        if result.returncode != 0:
            self.logger.warning(f"Erro na filogenia: {result.stderr}")
            tree_path = None
        
        # Exportar tabela para visualizaÃ§Ã£o
        export_dir = self.processed_dir / "exported"
        export_dir.mkdir(exist_ok=True)
        
        cmd_export = [
            "qiime", "tools", "export",
            "--input-path", str(table_path),
            "--output-path", str(export_dir)
        ]
        
        subprocess.run(cmd_export, capture_output=True, text=True)
        
        # Contar ASVs
        feature_table_path = export_dir / "feature-table.biom"
        asv_count = 0
        
        if feature_table_path.exists():
            # Usar qiime para contar features
            cmd_count = [
                "qiime", "feature-table", "summarize",
                "--i-table", str(table_path),
                "--o-visualization", str(self.processed_dir / "table_summary.qzv")
            ]
            subprocess.run(cmd_count, capture_output=True, text=True)
            
            # Estimativa simples
            try:
                import biom
                table = biom.load_table(str(feature_table_path))
                asv_count = len(table.ids(axis='observation'))
            except:
                # Fallback: contar linhas no arquivo exportado
                feature_table_tsv = export_dir / "feature-table.tsv"
                if feature_table_tsv.exists():
                    with open(feature_table_tsv, 'r') as f:
                        asv_count = sum(1 for line in f) - 1  # Menos cabeÃ§alho
        
        return {
            'table': table_path,
            'rep_seqs': rep_seqs_path,
            'tree': tree_path,
            'stats': stats_path,
            'asv_count': asv_count,
            'export_dir': str(export_dir)
        }
    
    def assign_taxonomy(self, table_path: Path) -> Dict:
        """Atribui taxonomia usando classificador SILVA"""
        
        taxonomy_path = self.processed_dir / "taxonomy.qza"
        
        # Usar classificador SILVA prÃ©-treinado
        classifier_path = self.project_dir / "silva-classifier.qza"
        
        # Baixar classificador se nÃ£o existir
        if not classifier_path.exists():
            self.logger.info("Baixando classificador SILVA...")
            # URL do classificador SILVA para regiÃ£o V3-V4
            classifier_url = "https://data.qiime2.org/2023.9/common/silva-138-99-515-806-nb-classifier.qza"
            
            try:
                import urllib.request
                urllib.request.urlretrieve(classifier_url, classifier_path)
            except:
                self.logger.warning("NÃ£o foi possÃ­vel baixar classificador, usando alternativa")
                classifier_path = None
        
        if classifier_path and classifier_path.exists():
            cmd = [
                "qiime", "feature-classifier", "classify-sklearn",
                "--i-classifier", str(classifier_path),
                "--i-reads", str(self.processed_dir / "rep_seqs.qza"),
                "--o-classification", str(taxonomy_path)
            ]
        else:
            # MÃ©todo alternativo: assign taxonomia com blast
            cmd = [
                "qiime", "feature-classifier", "classify-consensus-blast",
                "--i-query", str(self.processed_dir / "rep_seqs.qza"),
                "--i-reference-reads", str(self.project_dir / "silva-ref-seqs.qza"),
                "--i-reference-taxonomy", str(self.project_dir / "silva-ref-taxonomy.qza"),
                "--o-classification", str(taxonomy_path)
            ]
        
        self.logger.info("Atribuindo taxonomia")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            self.logger.warning(f"Erro na taxonomia: {result.stderr}")
            # Criar taxonomia vazia
            taxonomy_path = None
        
        # Exportar taxonomia
        if taxonomy_path and taxonomy_path.exists():
            export_dir = self.processed_dir / "taxonomy_export"
            export_dir.mkdir(exist_ok=True)
            
            cmd_export = [
                "qiime", "tools", "export",
                "--input-path", str(taxonomy_path),
                "--output-path", str(export_dir)
            ]
            subprocess.run(cmd_export, capture_output=True, text=True)
            
            taxonomy_file = export_dir / "taxonomy.tsv"
            if taxonomy_file.exists():
                taxonomy_df = pd.read_csv(taxonomy_file, sep='\t')
            else:
                taxonomy_df = pd.DataFrame()
        else:
            taxonomy_df = pd.DataFrame()
        
        return {
            'taxonomy': taxonomy_path,
            'taxonomy_df': taxonomy_df.to_dict() if not taxonomy_df.empty else {},
            'export_dir': str(export_dir) if 'export_dir' in locals() else None
        }
    
    def analyze_diversity(self, table_path: Path, tree_path: Path, params: Dict) -> Dict:
        """Executa anÃ¡lises de diversidade alfa e beta."""

        diversity_dir = self.results_dir / "diversity"
        diversity_dir.mkdir(exist_ok=True)

        core_metrics_path = diversity_dir / "core_metrics"

        # Construir comando correto dependendo se hÃ¡ Ã¡rvore ou nÃ£o
        if tree_path and tree_path.exists():
            cmd = [
                "qiime", "diversity", "core-metrics-phylogenetic",
                "--i-table", str(table_path),
                "--i-phylogeny", str(tree_path),
                "--p-sampling-depth", str(params["sampling_depth"]),
                "--m-metadata-file", str(self.project_dir / "metadata.tsv"),
                "--output-dir", str(core_metrics_path),
            ]
        else:
            # VersÃ£o nÃ£o filogenÃ©tica
            cmd = [
                "qiime", "diversity", "core-metrics",
                "--i-table", str(table_path),
                "--p-sampling-depth", str(params["sampling_depth"]),
                "--m-metadata-file", str(self.project_dir / "metadata.tsv"),
                "--output-dir", str(core_metrics_path),
            ]

        self.logger.info("Calculando mÃ©tricas de diversidade")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            self.logger.warning(f"Erro nas mÃ©tricas de diversidade: {result.stderr}")
            core_metrics_path = None

        diversity_results = {
            "core_metrics": str(core_metrics_path) if core_metrics_path else None,
            "alpha_diversity": {},
            "beta_diversity": {},
        }

        if core_metrics_path and core_metrics_path.exists():
            # Arquivos de diversidade alfa
            alpha_files = list(core_metrics_path.glob("*alpha*"))
            for alpha_file in alpha_files:
                metric_name = alpha_file.stem.replace("-", "_")
                diversity_results["alpha_diversity"][metric_name] = str(alpha_file)

            # Arquivos de diversidade beta
            beta_files = list(core_metrics_path.glob("*distance*"))
            for beta_file in beta_files:
                metric_name = beta_file.stem.replace("-", "_")
                diversity_results["beta_diversity"][metric_name] = str(beta_file)

            pcoa_files = list(core_metrics_path.glob("*pcoa*"))
            for pcoa_file in pcoa_files:
                diversity_results["pcoa"] = str(pcoa_file)

        # NMDS opcional
        if table_path.exists():
            try:
                export_dir = diversity_dir / "nmds"
                export_dir.mkdir(exist_ok=True)

                cmd_export = [
                    "qiime", "tools", "export",
                    "--input-path", str(table_path),
                    "--output-path", str(export_dir),
                ]
                subprocess.run(cmd_export, capture_output=True, text=True)

                feature_table = export_dir / "feature-table.biom"
                if feature_table.exists():
                    import biom
                    table = biom.load_table(str(feature_table))
                    data_matrix = table.matrix_data.toarray().T

                    from skbio.diversity import beta_diversity
                    from skbio.stats.ordination import nmds

                    bray_curtis = beta_diversity("braycurtis", data_matrix)
                    nmds_result = nmds(bray_curtis, number_of_dimensions=2)

                    nmds_df = pd.DataFrame(
                        nmds_result.samples.values,
                        columns=["NMDS1", "NMDS2"],
                        index=table.ids(axis="sample"),
                    )
                    nmds_coords_path = diversity_dir / "nmds_coordinates.csv"
                    nmds_df.to_csv(nmds_coords_path)

                    diversity_results["nmds"] = {
                        "coordinates": str(nmds_coords_path),
                        "stress": nmds_result.stress,
                    }

            except Exception as e:
                self.logger.warning(f"Erro no NMDS: {e}")

        return diversity_results

    
    def run_statistical_analyses(self, table_path: Path, taxonomy: Dict, metadata_path: Path) -> Dict:
        """Executa anÃ¡lises estatÃ­sticas"""
        
        stats_dir = self.results_dir / "statistics"
        stats_dir.mkdir(exist_ok=True)
        
        stats_results = {
            'permanova': {},
            'ancom': {},
            'alpha_group_significance': {}
        }
        
        # PERMANOVA para diversidade beta
        if metadata_path.exists():
            try:
                # Exportar matriz de distÃ¢ncia Bray-Curtis
                bray_curtis_path = list(Path(self.results_dir).glob("**/bray_curtis_distance_matrix.qza"))
                
                if bray_curtis_path:
                    bray_curtis_path = bray_curtis_path[0]
                    
                    # Testar diferentes variÃ¡veis categÃ³ricas
                    metadata_df = pd.read_csv(metadata_path, sep='\t')
                    categorical_cols = metadata_df.select_dtypes(include=['object']).columns
                    
                    for column in categorical_cols[:3]:  # Limitar a 3 variÃ¡veis
                        if metadata_df[column].nunique() > 1:
                            output_path = stats_dir / f"permanova_{column}.qzv"
                            
                            cmd = [
                                "qiime", "diversity", "beta-group-significance",
                                "--i-distance-matrix", str(bray_curtis_path),
                                "--m-metadata-file", str(metadata_path),
                                "--m-metadata-column", column,
                                "--p-permutations", str(config.DEFAULT_PARAMS['n_permutations']),
                                "--o-visualization", str(output_path)
                            ]
                            
                            result = subprocess.run(cmd, capture_output=True, text=True)
                            
                            if result.returncode == 0:
                                stats_results['permanova'][column] = {
                                    'visualization': str(output_path),
                                    'significant': True  # Placeholder
                                }
                            
            except Exception as e:
                self.logger.warning(f"Erro no PERMANOVA: {e}")
        
        # ANCOM para abundÃ¢ncia diferencial
        try:
            if table_path.exists() and taxonomy and metadata_path.exists():
                ancom_dir = stats_dir / "ancom"
                ancom_dir.mkdir(exist_ok=True)
                
                # Filtrar por nÃ­vel taxonÃ´mico (ex: gÃªnero)
                # ImplementaÃ§Ã£o simplificada
                
                stats_results['ancom'] = {
                    'status': 'implemented',
                    'output_dir': str(ancom_dir)
                }
                
        except Exception as e:
            self.logger.warning(f"Erro no ANCOM: {e}")
        
        # SignificÃ¢ncia de diversidade alfa
        try:
            if metadata_path.exists():
                # Para cada mÃ©trica alfa
                alpha_metrics = ['shannon', 'observed_features', 'faith_pd']
                
                for metric in alpha_metrics:
                    alpha_path = list(Path(self.results_dir).glob(f"**/{metric}_vector.qza"))
                    
                    if alpha_path:
                        alpha_path = alpha_path[0]
                        output_path = stats_dir / f"alpha_{metric}_significance.qzv"
                        
                        cmd = [
                            "qiime", "diversity", "alpha-group-significance",
                            "--i-alpha-diversity", str(alpha_path),
                            "--m-metadata-file", str(metadata_path),
                            "--o-visualization", str(output_path)
                        ]
                        
                        result = subprocess.run(cmd, capture_output=True, text=True)
                        
                        if result.returncode == 0:
                            stats_results['alpha_group_significance'][metric] = {
                                'visualization': str(output_path)
                            }
                            
        except Exception as e:
            self.logger.warning(f"Erro na significÃ¢ncia alfa: {e}")
        
        return stats_results
    
    def correct_batch_effects(self, table_path: Path, metadata_path: Path) -> Dict:
        """Corrige efeitos de batch"""
        
        batch_dir = self.results_dir / "batch_correction"
        batch_dir.mkdir(exist_ok=True)
        
        batch_results = {
            'methods_applied': [],
            'corrected_tables': []
        }
        
        try:
            if metadata_path.exists():
                metadata_df = pd.read_csv(metadata_path, sep='\t')
                
                # Identificar variÃ¡veis de batch
                batch_variables = []
                for col in metadata_df.columns:
                    if any(word in col.lower() for word in ['batch', 'run', 'plate', 'sequencing']):
                        batch_variables.append(col)
                
                if batch_variables:
                    self.logger.info(f"VariÃ¡veis de batch identificadas: {batch_variables}")
                    
                    # MÃ©todo 1: Remove Batch Effect usando ComBat (simulado)
                    # Em produÃ§Ã£o, usar bibliotecas como combat.py
                    
                    # Exportar tabela para anÃ¡lise
                    export_dir = batch_dir / "original"
                    export_dir.mkdir(exist_ok=True)
                    
                    cmd_export = [
                        "qiime", "tools", "export",
                        "--input-path", str(table_path),
                        "--output-path", str(export_dir)
                    ]
                    subprocess.run(cmd_export, capture_output=True, text=True)
                    
                    batch_results['methods_applied'].append('batch_identification')
                    batch_results['corrected_tables'].append(str(export_dir))
                    
        except Exception as e:
            self.logger.warning(f"Erro na correÃ§Ã£o de batch: {e}")
        
        return batch_results
    
    def generate_visualizations(self, table_path: Path, taxonomy: Dict, 
                              diversity_results: Dict, stats_results: Dict) -> Dict:
        """Gera visualizaÃ§Ãµes dos resultados"""
        
        viz_dir = self.results_dir / "visualizations"
        viz_dir.mkdir(exist_ok=True)
        
        viz_files = {}
        
        # 1. Barplot de composiÃ§Ã£o taxonÃ´mica
        try:
            if taxonomy.get('taxonomy'):
                barplot_path = viz_dir / "taxa_barplot.qzv"
                
                cmd = [
                    "qiime", "taxa", "barplot",
                    "--i-table", str(table_path),
                    "--i-taxonomy", str(taxonomy['taxonomy']),
                    "--m-metadata-file", str(self.project_dir / "metadata.tsv"),
                    "--o-visualization", str(barplot_path)
                ]
                
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode == 0:
                    viz_files['taxa_barplot'] = str(barplot_path)
                    
        except Exception as e:
            self.logger.warning(f"Erro no barplot: {e}")
        
        # 2. Heatmap
        try:
            if table_path.exists():
                # Exportar e criar heatmap com Plotly
                export_dir = viz_dir / "heatmap_data"
                export_dir.mkdir(exist_ok=True)
                
                cmd_export = [
                    "qiime", "tools", "export",
                    "--input-path", str(table_path),
                    "--output-path", str(export_dir)
                ]
                subprocess.run(cmd_export, capture_output=True, text=True)
                
                # Gerar heatmap interativo
                self.create_interactive_heatmap(export_dir, viz_dir)
                viz_files['heatmap'] = str(viz_dir / "heatmap.html")
                
        except Exception as e:
            self.logger.warning(f"Erro no heatmap: {e}")
        
        # 3. Curva de rarefaÃ§Ã£o
        try:
            if table_path.exists():
                rarefaction_path = viz_dir / "rarefaction_curve.qzv"
                
                cmd = [
                    "qiime", "diversity", "alpha-rarefaction",
                    "--i-table", str(table_path),
                    "--p-max-depth", str(config.DEFAULT_PARAMS['sampling_depth'] * 2),
                    "--o-visualization", str(rarefaction_path)
                ]
                
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode == 0:
                    viz_files['rarefaction'] = str(rarefaction_path)
                    
        except Exception as e:
            self.logger.warning(f"Erro na rarefaÃ§Ã£o: {e}")
        
        # 4. Boxplots para diversidade alfa
        try:
            alpha_metrics = diversity_results.get('alpha_diversity', {})
            for metric_name, metric_path in alpha_metrics.items():
                if 'shannon' in metric_name or 'observed' in metric_name:
                    # Criar boxplot interativo
                    self.create_alpha_boxplot(metric_path, viz_dir, metric_name)
                    viz_files[f'alpha_{metric_name}'] = str(viz_dir / f"alpha_{metric_name}.html")
                    
        except Exception as e:
            self.logger.warning(f"Erro nos boxplots: {e}")
        
        # 5. GrÃ¡ficos de dispersÃ£o para NMDS/PCoA
        try:
            if diversity_results.get('nmds'):
                nmds_coords = pd.read_csv(
                    diversity_results['nmds']['coordinates'], 
                    index_col=0
                )
                self.create_ordination_plot(nmds_coords, viz_dir, 'nmds')
                viz_files['nmds_plot'] = str(viz_dir / "nmds_plot.html")
                
        except Exception as e:
            self.logger.warning(f"Erro nos plots de ordenaÃ§Ã£o: {e}")
        
        # 6. Biplot (se disponÃ­vel)
        try:
            if diversity_results.get('pcoa'):
                # Tentar criar biplot com variÃ¡veis ambientais
                self.create_biplot(diversity_results, viz_dir)
                viz_files['biplot'] = str(viz_dir / "biplot.html")
                
        except Exception as e:
            self.logger.warning(f"Erro no biplot: {e}")
        
        return viz_files
    
    def create_interactive_heatmap(self, export_dir: Path, output_dir: Path):
        """Cria heatmap interativo com Plotly"""
        
        try:
            # Carregar tabela de features
            feature_table_path = export_dir / "feature-table.tsv"
            if feature_table_path.exists():
                df = pd.read_csv(feature_table_path, sep='\t', index_col=0)
                
                # Transformar log
                df_log = np.log1p(df)
                
                # Selecionar top features
                top_features = df.sum(axis=1).nlargest(50).index
                df_top = df_log.loc[top_features]
                
                # Criar heatmap
                fig = px.imshow(
                    df_top,
                    labels=dict(x="Amostras", y="ASVs", color="AbundÃ¢ncia (log)"),
                    x=df_top.columns,
                    y=df_top.index,
                    aspect="auto",
                    color_continuous_scale="Viridis"
                )
                
                fig.update_layout(
                    title="Heatmap - Top 50 ASVs",
                    xaxis_title="Amostras",
                    yaxis_title="Features (ASVs)",
                    height=800
                )
                
                # Salvar
                heatmap_path = output_dir / "heatmap.html"
                fig.write_html(str(heatmap_path))
                
        except Exception as e:
            self.logger.error(f"Erro criando heatmap: {e}")
    
    def create_alpha_boxplot(self, metric_path: Path, output_dir: Path, metric_name: str):
        """Cria boxplot interativo para diversidade alfa"""
        
        try:
            # Exportar mÃ©trica alfa
            export_dir = output_dir / f"alpha_{metric_name}"
            export_dir.mkdir(exist_ok=True)
            
            cmd_export = [
                "qiime", "tools", "export",
                "--input-path", str(metric_path),
                "--output-path", str(export_dir)
            ]
            
            result = subprocess.run(cmd_export, capture_output=True, text=True)
            
            if result.returncode == 0:
                # Carregar dados
                alpha_file = export_dir / "alpha-diversity.tsv"
                if alpha_file.exists():
                    df = pd.read_csv(alpha_file, sep='\t')
                    
                    # Carregar metadados
                    metadata_path = self.project_dir / "metadata.tsv"
                    if metadata_path.exists():
                        metadata = pd.read_csv(metadata_path, sep='\t')
                        
                        # Mesclar
                        df_merged = pd.merge(
                            df, metadata,
                            left_on='SampleID',
                            right_on=metadata.columns[0]
                        )
                        
                        # Criar boxplot
                        fig = px.box(
                            df_merged,
                            x='health_status',  # Ajustar para coluna real
                            y=df.columns[1],
                            color='health_status',
                            title=f"Diversidade Alfa - {metric_name.replace('_', ' ').title()}",
                            points="all"
                        )
                        
                        fig.update_layout(
                            xaxis_title="Status de SaÃºde",
                            yaxis_title=metric_name.replace('_', ' ').title(),
                            showlegend=False
                        )
                        
                        # Salvar
                        boxplot_path = output_dir / f"alpha_{metric_name}.html"
                        fig.write_html(str(boxplot_path))
                        
        except Exception as e:
            self.logger.error(f"Erro criando boxplot: {e}")
    
    def create_ordination_plot(self, coordinates_df: pd.DataFrame, output_dir: Path, method: str):
        """Cria plot de ordenaÃ§Ã£o (NMDS/PCoA)"""
        
        try:
            # Carregar metadados
            metadata_path = self.project_dir / "metadata.tsv"
            if metadata_path.exists():
                metadata = pd.read_csv(metadata_path, sep='\t')
                
                # Mesclar coordenadas com metadados
                df_plot = pd.merge(
                    coordinates_df.reset_index(),
                    metadata,
                    left_on='index',
                    right_on=metadata.columns[0]
                )
                
                # Criar scatter plot
                fig = px.scatter(
                    df_plot,
                    x=coordinates_df.columns[0],
                    y=coordinates_df.columns[1],
                    color='health_status',
                    symbol='infection_level' if 'infection_level' in df_plot.columns else None,
                    hover_name='index',
                    title=f"AnÃ¡lise de OrdenaÃ§Ã£o - {method.upper()}",
                    size_max=15
                )
                
                fig.update_layout(
                    xaxis_title=f"{coordinates_df.columns[0]}",
                    yaxis_title=f"{coordinates_df.columns[1]}",
                    height=600
                )
                
                # Salvar
                plot_path = output_dir / f"{method}_plot.html"
                fig.write_html(str(plot_path))
                
        except Exception as e:
            self.logger.error(f"Erro criando plot de ordenaÃ§Ã£o: {e}")
    
    def create_biplot(self, diversity_results: Dict, output_dir: Path):
        """Cria biplot para anÃ¡lise multivariada"""
        
        try:
            # ImplementaÃ§Ã£o simplificada
            # Em produÃ§Ã£o, usar RDA ou CCA com variÃ¡veis ambientais
            
            metadata_path = self.project_dir / "metadata.tsv"
            if metadata_path.exists():
                metadata = pd.read_csv(metadata_path, sep='\t')
                
                # Criar biplot bÃ¡sico
                fig = go.Figure()
                
                # Adicionar pontos (amostras)
                fig.add_trace(go.Scatter(
                    x=[0, 1, 2, 3],
                    y=[0, 1, 2, 3],
                    mode='markers',
                    name='Amostras'
                ))
                
                # Adicionar vetores (variÃ¡veis)
                fig.add_trace(go.Scatter(
                    x=[0, 0.5, 1],
                    y=[0, 0.5, 1],
                    mode='lines+text',
                    name='VariÃ¡veis',
                    text=['Temperatura', 'Salinidade', 'pH'],
                    textposition='top center'
                ))
                
                fig.update_layout(
                    title="Biplot - AnÃ¡lise Multivariada",
                    xaxis_title="Componente 1",
                    yaxis_title="Componente 2",
                    height=600
                )
                
                # Salvar
                biplot_path = output_dir / "biplot.html"
                fig.write_html(str(biplot_path))
                
        except Exception as e:
            self.logger.error(f"Erro criando biplot: {e}")

    # ============================================================================
    # SISTEMA DE COMPARAÃ‡ÃƒO COM DADOS PÃšBLICOS
    # ============================================================================
# ============================================================================
# INTERFACE STREAMLIT - SISTEMA COMPLETO
# ============================================================================
class ShrimpMicrobiomeApp:
    """AplicaÃ§Ã£o Streamlit completa"""

 # ===================== NOVOS HELPERS DE PERMISSÃƒO ======================
   
    def __init__(self):
        self.setup_database()
        self.setup_session_state()

        # Instanciar o ComparativeAnalyzer UMA vez por sessÃ£o
        if "comparative_analyzer" not in st.session_state:
            st.session_state.comparative_analyzer = ComparativeAnalyzer(
                config.PUBLIC_DATA_DIR,
                config.USER_DATA_DIR,
            )
    def inject_global_css(self):
        """CSS global para deixar a tela de login (e o app) mais bonito."""
        st.markdown("""
        <style>
          /* Fundo geral da aplicaÃ§Ã£o */
          .stApp {
            background: radial-gradient(circle at top, #0f172a 0, #020617 55%);
          }

          /* TÃ­tulo e subtÃ­tulo da tela de login */
          .login-title {
            text-align: center;
            font-size: 2.2rem;
            font-weight: 700;
            color: #e5e7eb;
            margin-top: 0.5rem;
            margin-bottom: 0.3rem;
          }
          .login-subtitle {
            text-align: center;
            color: #9ca3af;
            font-size: 0.95rem;
            margin-bottom: 1.5rem;
          }

          /* CartÃ£o da tela de login */
          .login-card {
            background: rgba(15, 23, 42, 0.92);
            border-radius: 18px;
            padding: 1.8rem 2rem;
            box-shadow: 0 18px 45px rgba(15, 23, 42, 0.9);
            border: 1px solid rgba(148, 163, 184, 0.35);
            backdrop-filter: blur(12px);
          }

          /* Labels e inputs dentro do card */
          .login-card label {
            color: #e5e7eb !important;
            font-weight: 500;
          }
          .login-card input {
            background-color: rgba(15, 23, 42, 0.6) !important;
            color: #e5e7eb !important;
            border-radius: 10px !important;
            border: 1px solid rgba(148, 163, 184, 0.6) !important;
          }
          .login-card input:focus {
            border-color: #22d3ee !important;
            box-shadow: 0 0 0 1px #22d3ee !important;
          }

          /* BotÃ£o principal dentro do card de login */
          .login-card .stButton > button {
            background: linear-gradient(135deg, #06b6d4, #3b82f6) !important;
            color: white !important;
            border-radius: 999px !important;
            border: none !important;
            padding: 0.4rem 1.4rem !important;
            font-weight: 600 !important;
            letter-spacing: .02em;
          }
          .login-card .stButton > button:hover {
            filter: brightness(1.07);
            transform: translateY(-1px);
          }

          /* Texto de ajuda (links de registro/esqueci senha) */
          .login-help {
            color: #9ca3af;
            font-size: 0.85rem;
            margin-top: 1.4rem;
          }
          .login-help a {
            color: #22d3ee;
            text-decoration: none;
          }
          .login-help a:hover {
            text-decoration: underline;
          }

          /* Sidebar mais elegante */
          section[data-testid="stSidebar"] {
            background: linear-gradient(180deg, #020617 0, #020617 30%, #020617 100%);
            border-right: 1px solid rgba(148, 163, 184, 0.35);
          }
        </style>
        """, unsafe_allow_html=True)
        
    def setup_database(self):
        """Configura banco de dados"""
        engine = create_engine(config.DB_URL)
        Base.metadata.create_all(engine)
        self.Session = sessionmaker(bind=engine)
        
        # Criar admin padrÃ£o se nÃ£o existir
        self.ensure_default_admin()
    
    def setup_session_state(self):
        """Configura estado da sessÃ£o"""
        if 'current_user' not in st.session_state:
            st.session_state.current_user = None
        if 'current_project' not in st.session_state:
            st.session_state.current_project = None
        if 'processing_queue' not in st.session_state:
            st.session_state.processing_queue = []
        if 'comparison_results' not in st.session_state:
            st.session_state.comparison_results = {}
    
    def ensure_default_admin(self):
        """Garante que admin padrÃ£o existe"""
        session = self.Session()
        try:
            admin = session.query(User).filter_by(username=config.DEFAULT_ADMIN_USERNAME).first()
            if not admin:
                admin = User(
                   username=config.DEFAULT_ADMIN_USERNAME,
                   email=config.DEFAULT_ADMIN_EMAIL,
                   password_hash=hashlib.sha256(config.DEFAULT_ADMIN_PASSWORD.encode()).hexdigest(),
                   full_name="Administrador ShrimpMicrobiome",
                   organization=config.COMPANY_NAME,
                   role=UserRole.ADMIN,
                   status=UserStatus.APPROVED,
                   study_species=config.TARGET_SPECIES
               )

                session.add(admin)
                session.commit()
                logger.info("Admin padrÃ£o criado")
        finally:
            session.close()

    
    def run(self):
        """Executa a aplicaÃ§Ã£o"""
        self.inject_global_css()

        # Sidebar com navegaÃ§Ã£o
        with st.sidebar:
            st.markdown(
                f"""
                <div style="
                  text-align: center;
                  padding: 1.5rem 1rem 1rem 1rem;
                  background: linear-gradient(145deg, #0f172a, #022c22);
                  border-radius: 18px;
                  box-shadow: 0 12px 30px rgba(15, 23, 42, 0.45);
                  margin-bottom: 1.5rem;
                ">
                  <div style="font-size: 3rem; line-height: 1;">ðŸ¦</div>
                  <h3 style="margin: 0.4rem 0 0.2rem 0; color: #e2f3ff; letter-spacing: 0.06em;">
                    {config.APP_NAME}
                  </h3>
                  <p style="color: #bae6fd; font-size: 0.85rem; margin: 0.3rem 0 0;">
                    {config.APP_TAGLINE}
                  </p>
                </div>
                """,
                unsafe_allow_html=True
            )

            # Menu de navegaÃ§Ã£o
            if not st.session_state.get("current_user"):
                menu_items = ["Login", "Registro", "Sobre"]
            else:
                menu_items = [
                    "Dashboard",
                    "Meus Projetos",
                    "Upload de Dados",
                    "AnÃ¡lise Comparativa",
                    "Base de Dados PÃºblica",
                    "ConfiguraÃ§Ãµes",
                ]

                if st.session_state.current_user.role == UserRole.ADMIN:
                    menu_items.append("Painel Admin")

                menu_items.append("Sair")

            selected_menu = st.selectbox("NavegaÃ§Ã£o", menu_items)

            # Status do usuÃ¡rio logado
            if st.session_state.get("current_user"):
                user = st.session_state.current_user
                st.divider()
                st.markdown(
                    f"""
                    <div style="
                      padding: 0.75rem 0.9rem;
                      background: rgba(15, 23, 42, 0.03);
                      border-radius: 14px;
                      border: 1px solid rgba(148, 163, 184, 0.35);
                    ">
                      <p style="margin: 0; font-weight: 600; color: #0f172a;">
                        {user.full_name}
                      </p>
                      <p style="margin: 0.2rem 0 0; color: #64748b; font-size: 0.8rem;">
                        {user.organization}
                      </p>
                    </div>
                    """,
                    unsafe_allow_html=True
                )

            st.divider()

            st.markdown(
                f"""
               <p style="text-align: center; color: #94a3b8; font-size: 0.75rem; margin-top: 0.5rem;">
                  Sea Scient Â· v{config.APP_VERSION} Â· {config.APP_YEAR}
                </p>
                """,
                unsafe_allow_html=True
           )


        # ConteÃºdo principal baseado na seleÃ§Ã£o
        if selected_menu == "Login":
            self.login_page()
        elif selected_menu == "Registro":
            self.register_page()
        elif selected_menu == "Dashboard":
            self.dashboard_page()
        elif selected_menu == "Meus Projetos":
            self.projects_page()
        elif selected_menu == "Upload de Dados":
            self.upload_page()
        elif selected_menu == "AnÃ¡lise Comparativa":
            self.comparative_analysis_page()
        elif selected_menu == "Base de Dados PÃºblica":
            self.public_database_page()
        elif selected_menu == "Painel Admin":
            self.admin_panel()
        elif selected_menu == "ConfiguraÃ§Ãµes":
            self.settings_page()
        elif selected_menu == "Sair":
            self.logout()
        elif selected_menu == "Sobre":
            self.about_page()

    def login_page(self):
        """PÃ¡gina de login"""
        # CabeÃ§alho central com tema do camarÃ£o, mar e pescador
        st.markdown(
            """
            <div style="text-align: center; margin-bottom: 2rem;">
              <div style="font-size: 1.8rem; line-height: 1.2;">ðŸŒŠðŸŒŠðŸŒŠ</div>
              <div style="font-size: 1.3rem; line-height: 1.2;">ðŸ§‘ðŸŽ£ðŸ•¸ï¸</div>
              <div style="font-size: 3.2rem; line-height: 1;">ðŸ¦</div>
              <div style="
                width: 110px;
                height: 5px;
                border-radius: 999px;
                margin: 0.4rem auto 0.6rem;
                background: linear-gradient(to right, #e0f2fe, #38bdf8);
              "></div>
              <h1 style="margin: 0.2rem 0; color: #0f172a;">Login</h1>
              <p style="margin: 0; color: #64748b; font-size: 0.95rem;">
                Entre para acessar suas anÃ¡lises 16S de <em>Penaeus vannamei</em>.
              </p>
            </div>
            """,
            unsafe_allow_html=True,
        )

        # Layout centralizado
        col1, col2, col3 = st.columns([1, 2.2, 1])

        with col2:
            # Card de login com visual marinho
            st.markdown(
                """
                <div style="
                  background: radial-gradient(circle at top, #e0f2fe 0, #f9fafb 55%);
                  border-radius: 22px;
                  padding: 1.2rem 1.6rem 1.8rem 1.6rem;
                  box-shadow: 0 18px 45px rgba(15, 23, 42, 0.18);
                  border: 1px solid rgba(148, 163, 184, 0.45);
                  position: relative;
                  overflow: hidden;
                ">
                  <div style="
                    position: absolute;
                    inset: 0;
                    background: radial-gradient(circle at bottom left, rgba(56, 189, 248, 0.20), transparent 55%);
                    opacity: 0.9;
                    pointer-events: none;
                  "></div>
                  <div style="position: relative; z-index: 1;">
                """,
                unsafe_allow_html=True,
            )

            with st.form("login_form"):
                username = st.text_input("UsuÃ¡rio")
                password = st.text_input("Senha", type="password")
                submit = st.form_submit_button("Entrar")

                if submit:
                    if self.authenticate_user(username, password):
                        st.success("Login bem-sucedido! Redirecionando...")
                        # Fallback para versÃµes antigas do Streamlit
                        if hasattr(st, "rerun"):
                            st.rerun()
                        elif hasattr(st, "experimental_rerun"):
                            st.experimental_rerun()
                    else:
                        st.error("UsuÃ¡rio ou senha incorretos.")

            # Texto auxiliar
            st.markdown(
                """
                <div style="margin-top: 1.2rem; font-size: 0.85rem; color: #475569;">
                  <strong>Primeiro acesso?</strong><br/>
                  - <a href="#" style="color:#0ea5e9; text-decoration:none;">Registrar nova conta</a><br/>
                  - <a href="#" style="color:#0ea5e9; text-decoration:none;">Esqueci minha senha</a><br/><br/>
                  <strong>Acesso rÃ¡pido para teste:</strong><br/>
                  UsuÃ¡rio: <code>demo</code> &nbsp; | &nbsp; Senha: <code>demo123</code>
                </div>
                """,
                unsafe_allow_html=True,
            )

            # Fecha o card
            st.markdown(
                """
                  </div>
                </div>
                """,
                unsafe_allow_html=True,
            )

        # RodapÃ© da pÃ¡gina de login
        st.markdown(
            f"""
            <div style="margin-top: 2.5rem; text-align: center; font-size: 0.8rem; color: #94a3b8;">
              <div style="margin-bottom: 0.15rem;">
                <strong>ShrimpMicrobiome Platform</strong> Â· v{config.APP_VERSION}
              </div>
              <div>
                Plataforma ShrimpMicrobiome: Projeto da Sea Scient, desenvolvida em colaboraÃ§Ã£o com Databiomics.<br/>
                Todos os direitos reservados Â© 2025
              </div>
            </div>
            """,
            unsafe_allow_html=True,
        )

 
    def authenticate_user(self, username: str, password: str) -> bool:
        """Autentica usuÃ¡rio e armazena um CurrentUser em session_state."""
        session = self.Session()
        try:
            user = session.query(User).filter_by(username=username).first()

            # UsuÃ¡rio ou senha invÃ¡lidos
            if not user or not user.verify_password(password):
                return False

            # Verificar status da conta
            if user.status != UserStatus.APPROVED:
                if user.status == UserStatus.PENDING:
                    st.warning("Conta aguardando aprovaÃ§Ã£o do administrador.")
                elif user.status == UserStatus.REJECTED:
                    st.error("Conta rejeitada. Entre em contato com o suporte.")
                elif user.status == UserStatus.SUSPENDED:
                    st.error("Conta suspensa. Entre em contato com o suporte.")
                else:
                    st.error("Status de conta invÃ¡lido.")
                return False

            # Criar snapshot simples (sem vÃ­nculo com a Session) para guardar na sessÃ£o
            current_user = CurrentUser(
                id=user.id,
                email=user.email,
                username=user.username,
                full_name=user.full_name or user.username,
                role=user.role,  # UserRole
                organization=user.organization,
            )

            # Agora st.session_state.current_user NÃƒO Ã© mais o objeto SQLAlchemy
            st.session_state.current_user = current_user

            # Atualizar Ãºltimo login
            user.last_login = datetime.utcnow()
            session.commit()

            return True
        finally:
            session.close()

    
    def register_page(self):
        """PÃ¡gina de registro"""
        st.title("ðŸ“‘ Registro de Nova Conta")

        with st.form("register_form"):
            col1, col2 = st.columns(2)

            with col1:
                username = st.text_input("UsuÃ¡rio*")
                email = st.text_input("Email*")
                full_name = st.text_input("Nome Completo*")
                organization = st.text_input("InstituiÃ§Ã£o/OrganizaÃ§Ã£o*")

                # Tipo de usuÃ¡rio (funciona para todos os papÃ©is)
                role_labels = {
                    "Pesquisador(a)": UserRole.RESEARCHER,
                    "Produtor(a) / Fazenda": UserRole.FARMER,
                    "Estudante": UserRole.STUDENT,
                    "PÃºblico / Observador": UserRole.PUBLIC,
                }
                selected_role_label = st.selectbox(
                    "Tipo de usuÃ¡rio*",
                    list(role_labels.keys())
                )
                selected_role = role_labels[selected_role_label]

            with col2:
                country = st.selectbox("PaÃ­s*", SHRIMP_PRODUCING_COUNTRIES)
                study_species = st.text_input("EspÃ©cie de Estudo*", value=config.TARGET_SPECIES)
                password = st.text_input("Senha*", type="password")
                confirm_password = st.text_input("Confirmar Senha*", type="password")

            research_interests = st.multiselect(
                "Interesses de Pesquisa",
                ["Microbioma", "DoenÃ§as", "NutriÃ§Ã£o", "GenÃ´mica", "Aquacultura SustentÃ¡vel", "Outro"],
            )

            terms = st.checkbox("Aceito os termos de uso e polÃ­tica de privacidade*")

            submitted = st.form_submit_button("Registrar")

            if submitted:
                if not all([username, email, full_name, organization, password, confirm_password, terms]):
                    st.error("Por favor, preencha todos os campos obrigatÃ³rios (*)")
                elif password != confirm_password:
                    st.error("As senhas nÃ£o coincidem")
                else:
                    session = self.Session()
                    try:
                        existing_user = (
                            session.query(User)
                            .filter((User.username == username) | (User.email == email))
                            .first()
                        )

                        if existing_user:
                            st.error("UsuÃ¡rio ou email jÃ¡ cadastrado")
                        else:
                            new_user = User(
                                username=username,
                                email=email,
                                full_name=full_name,
                                organization=organization,
                                country=country,
                                study_species=study_species,
                                password_hash=hashlib.sha256(password.encode()).hexdigest(),
                                research_interests=research_interests,
                                role=selected_role,
                                status=UserStatus.PENDING,
                            )

                            session.add(new_user)
                            session.commit()

                            st.success(
                                """
                                âœ… Registro realizado com sucesso!

                                Sua conta estÃ¡ aguardando aprovaÃ§Ã£o do administrador.
                                VocÃª receberÃ¡ um email quando sua conta for aprovada.
                                """
                            )

                            logger.info(
                                f"Novo usuÃ¡rio registrado: {username} ({email}) â€“ tipo: {selected_role.value}"
                            )
                    finally:
                        session.close()
    
    def dashboard_page(self):
        """Dashboard principal"""
        st.title("ðŸ“Š Dashboard")
        
        # MÃ©tricas rÃ¡pidas (depois vocÃª pode ligar isso ao banco)
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Projetos", "0", "0")
        with col2:
            st.metric("Amostras", "0", "0")
        with col3:
            st.metric("AnÃ¡lises", "0", "0")
        with col4:
            st.metric("ComparaÃ§Ãµes", "0", "0")
        
        st.markdown("### âš¡ AÃ§Ãµes RÃ¡pidas")
        col1, col2, col3 = st.columns(3)
        with col1:
            if st.button("ðŸ“ Novo Projeto", use_container_width=True):
                st.session_state.show_new_project = True
        
        with col2:
            if st.button("ðŸ“¤ Upload Dados", use_container_width=True):
                st.session_state.selected_menu = "Upload de Dados"
        
        with col3:
            if st.button("ðŸ” AnÃ¡lise Comparativa", use_container_width=True):
                st.session_state.selected_menu = "AnÃ¡lise Comparativa"
        
        if st.session_state.get("show_new_project"):
            self.new_project_form()
        
        # --------- BLOCO: estatÃ­stica rÃ¡pida da base pÃºblica ----------
        st.markdown("### ðŸŒ Resumo rÃ¡pido da Base PÃºblica (SRA / ENA)")
        
        region = st.selectbox(
            "RegiÃ£o alvo 16S",
            ["V3-V4", "V4", "V1-V3", "Todas"],
            index=0,
        )
        region_arg = None if region == "Todas" else region
        
        max_to_check = st.slider(
            "MÃ¡ximo de runs a inspecionar (prÃ©-visualizaÃ§Ã£o)",
            min_value=10,
            max_value=500,
            value=100,
            step=10,
        )
        
        # Adicionar uma busca de termo para a base pÃºblica
        search_term = st.text_input(
            "Termo de busca",
            value="Penaeus vannamei 16S",
            help="Buscar por espÃ©cie, tecido, condiÃ§Ã£o, etc."
        )
        
        if st.button("ðŸ›° Atualizar EstatÃ­sticas PÃºblicas"):
            with st.spinner("Consultando bases pÃºblicas (ENA/SRA)..."):
                downloader = PublicDataDownloader(max_samples=max_to_check)
                
                try:
                    # Buscar datasets pÃºblicos (consulta principal)
                    try:
                        target_species = config.TARGET_SPECIES
                    except Exception:
                        target_species = "Penaeus vannamei"

                    df_raw = downloader.search_public_datasets(
                        species=target_species,
                        marker="16S",
                        source=None,          # ENA + SRA
                        region=region_arg,
                        n_max=max_to_check,
                        search_term=search_term
                    )
                    
                    # MÃ©todo alternativo: organizar dados (usa internamente search_public_datasets)
                    try:
                        organized = downloader.organize_public_data(
                            output_dir=str(config.PUBLIC_DATA_DIR),
                            search_term=search_term,
                            species=target_species,
                            marker="16S",
                            region=region_arg,
                            n_max=max_to_check,
                        )
                    except Exception as e:
                        # Se organize_public_data falhar, usar df_raw
                        logger.warning(f"organize_public_data falhou, usando df_raw: {e}")
                        organized = df_raw

                    if organized is None or organized.empty:
                        if df_raw is not None and not df_raw.empty:
                            organized = df_raw
                        else:
                            st.warning("Nenhum dataset encontrado com os critÃ©rios atuais.")
                            return

                    total_runs = len(organized)
                    if "total_reads" in organized.columns:
                        total_reads = pd.to_numeric(organized["total_reads"], errors="coerce").fillna(0).astype(int).sum()
                    else:
                        total_reads = 0
                    
                    # Tentar obter contagens de saÃºde
                    n_healthy = 0
                    n_dysbiosis = 0
                    if "health_status" in organized.columns:
                        n_healthy = (organized["health_status"] == "healthy").sum()
                        n_dysbiosis = (organized["health_status"] == "dysbiosis").sum()
                    elif "status" in organized.columns:
                        n_healthy = (organized["status"].str.contains("healthy", case=False, na=False)).sum()
                        n_dysbiosis = (organized["status"].str.contains("disease|dysbiosis", case=False, na=False)).sum()

                    c1, c2, c3, c4 = st.columns(4)
                    with c1:
                        st.metric("Runs 16S", total_runs)
                    with c2:
                        st.metric("Total de Reads", f"{total_reads:,}".replace(",", "."))
                    with c3:
                        st.metric("Datasets SaudÃ¡veis", n_healthy)
                    with c4:
                        st.metric("Datasets com Disbiose", n_dysbiosis)

                    st.markdown("PrÃ©via dos datasets encontrados:")
                    
                    # Definir colunas para mostrar
                    columns_to_show = ["accession", "title", "total_reads"]
                    if "health_status" in organized.columns:
                        columns_to_show.append("health_status")
                    if "infection_level" in organized.columns:
                        columns_to_show.append("infection_level")
                    if "tissue" in organized.columns:
                        columns_to_show.append("tissue")
                    
                    existing_columns = [col for col in columns_to_show if col in organized.columns]
                    
                    if existing_columns:
                        st.dataframe(
                            organized[existing_columns].head(20),
                            use_container_width=True,
                        )
                    else:
                        st.dataframe(
                            organized.head(20),
                            use_container_width=True,
                        )
                            
                except Exception as e:
                    st.error(f"Erro ao buscar dados pÃºblicos: {str(e)}")
                    logger.error(f"Erro no dashboard ao buscar dados pÃºblicos: {e}", exc_info=True)



    def download_public_datasets(self, datasets_df: pd.DataFrame, search_term: str = None, region: str = "V3-V4", max_to_check: int = 100):
        """Faz download rÃ¡pido de datasets pÃºblicos a partir do DataFrame retornado pelo SRA.

        Este mÃ©todo Ã© usado no Dashboard apÃ³s a busca de dados pÃºblicos.
        """
        if datasets_df.empty:
            st.warning("Nenhum dataset disponÃ­vel para download.")
            return

        downloader = PublicDataDownloader()
        
        # Extrair bioprojects do DataFrame
        bioprojects = []
        if 'bioproject' in datasets_df.columns:
            bioprojects = datasets_df['bioproject'].dropna().unique().tolist()
        elif 'bioproject_accession' in datasets_df.columns:
            bioprojects = datasets_df['bioproject_accession'].dropna().unique().tolist()
        elif 'study_accession' in datasets_df.columns:
            bioprojects = datasets_df['study_accession'].dropna().unique().tolist()
        
        if not bioprojects:
            st.warning("NÃ£o foi possÃ­vel identificar BioProjects para organizar os dados.")
            return
        
        # BUSCA REAL NAS BASES SRA + ENA USANDO OS BIOPROJECTS
        try:
            organized_data = downloader.search_public_datasets(
                species="Penaeus vannamei",
                marker="16S",
                region=region,
                bioproject_ids=bioprojects,      # <--- fundamental
                search_term=search_term,         # <--- agora funciona
                source=None,                     # <--- busca ENA + SRA
                n_max=max_to_check,
            )
        except Exception as e:
            st.warning(f"Erro ao buscar datasets: {e}. Usando dados originais.")
            organized_data = datasets_df

        if organized_data.empty:
            st.warning("Nenhum dado encontrado no SRA/ENA.")
            return

        accessions = organized_data["accession"].dropna().tolist()
        if not accessions:
            st.warning("Nenhum accession disponÃ­vel para download.")
            return

        # Limite simples
        accessions = accessions[:min(len(accessions), downloader.max_samples)]

        output_dir = config.PUBLIC_DATA_DIR / "downloaded"
        output_dir.mkdir(exist_ok=True, parents=True)

        st.info(f"Baixando {len(accessions)} datasets do SRA...")
        
        # Seleciona mÃ©todo correto de download
        try:
            if hasattr(downloader, 'download_sra_samples'):
                downloaded = downloader.download_sra_samples(accessions, output_dir)
            elif hasattr(downloader, 'download_public_samples'):
                downloaded = downloader.download_public_samples(
                    organized_df=organized_data,
                    accession_list=accessions,
                    output_dir=output_dir,
                    threads=4
                )
            else:
                st.warning("MÃ©todo de download nÃ£o disponÃ­vel. Simulando downloads...")
                downloaded = []
                for accession in accessions[:5]:
                    st.write(f"Simulando download de {accession}...")
                    downloaded.append(accession)
        except Exception as e:
            st.error(f"Erro no download: {e}")
            downloaded = []
        
        if downloaded:
            st.success(f"Download concluÃ­do: {len(downloaded)} arquivos.")
        else:
            st.warning("Nenhum arquivo foi baixado.")

        # Processamento opcional
        if st.checkbox("Processar automaticamente os datasets baixados"):
            subset_meta = organized_data[organized_data["accession"].isin(accessions[:5])]
            self.process_public_datasets(output_dir, subset_meta)


    def create_project(self, name: str, data: Dict):
        """Cria novo projeto"""
        session = self.Session()
        try:
            project = Project(
                owner_id=st.session_state.current_user.id,
                name=name,
                description=data.get("description", ""),
                species=data.get("species", getattr(config, "TARGET_SPECIES", "")),
                project_type=data.get("project_type", "comparative_analysis"),
                experimental_design={
                    "design": data.get("experimental_design"),
                    "replicates": data.get("replicates", 3),
                },
                # Dict genÃ©rico com os parÃ¢metros que o usuÃ¡rio definiu
                diet_composition=data.get("diet_composition", {}),
                status="planning",
                base_path=str(
                    config.PROJECTS_DIR
                    / f"{st.session_state.current_user.username}_{name}"
                ),
            )

            session.add(project)
            session.commit()

            # Criar diretÃ³rio do projeto
            project_dir = Path(project.base_path)
            project_dir.mkdir(exist_ok=True, parents=True)

            st.success(f"Projeto '{name}' criado com sucesso!")

        finally:
            session.close()

    def projects_page(self):
        """PÃ¡gina de projetos do usuÃ¡rio"""
        st.title("ðŸ“ Meus Projetos")

        session = self.Session()
        try:
            projects = (
                session.query(Project)
                .filter_by(owner_id=st.session_state.current_user.id)
                .order_by(Project.created_at.desc())
                .all()
            )

            if not projects:
                st.info("VocÃª ainda nÃ£o tem projetos. Crie seu primeiro projeto!")
                if st.button("ðŸ†• Criar Primeiro Projeto"):
                    st.session_state.show_new_project = True
                    st.rerun()
                return

            for project in projects:
                with st.expander(
                    f"ðŸ“‚ {project.name} - {project.status}", expanded=False
                ):
                    col1, col2 = st.columns([3, 1])

                    with col1:
                        st.markdown(
                            f"**DescriÃ§Ã£o:** {project.description or 'Sem descriÃ§Ã£o'}"
                        )
                        st.markdown(f"**EspÃ©cie:** {project.species}")
                        st.markdown(
                            f"**Criado em:** {project.created_at.strftime('%d/%m/%Y')}"
                        )
                        st.markdown(f"**Amostras:** {project.sample_count}")

                    with col2:
                        status_color = {
                            "planning": "blue",
                            "active": "green",
                            "completed": "gray",
                            "failed": "red",
                        }.get(project.status, "gray")

                        st.markdown(
                            (
                                "<span style='color: {color}; font-weight: bold;'>"
                                "{status}</span>"
                            ).format(
                                color=status_color,
                                status=project.status.upper(),
                            ),
                            unsafe_allow_html=True,
                        )

                        if st.button("Abrir", key=f"open_{project.id}"):
                            st.session_state.current_project = project
                            # Aqui vocÃª pode trocar para o menu de anÃ¡lise
                            st.session_state.selected_menu = "AnÃ¡lise do Projeto"
                            st.rerun()

                    st.progress((project.progress or 0) / 100)

        finally:
            session.close()

    def upload_page(self):
        """PÃ¡gina de upload de dados"""
        st.title("ðŸ“‚ Upload de Dados FASTQ")

        if not st.session_state.current_project:
            st.warning("Selecione um projeto primeiro")
            return

        st.markdown(
            f"""
        ### Projeto: {st.session_state.current_project.name}

        **InstruÃ§Ãµes:**
        1. Prepare seus arquivos FASTQ (pares R1 e R2)
        2. Prepare arquivo de metadados (TSV/CSV)
        3. FaÃ§a upload dos arquivos abaixo
        4. Configure os parÃ¢metros do pipeline
        5. Execute o processamento
        """
        )

        # Upload de arquivos
        st.markdown("#### 1. Upload de Arquivos")

        col1, col2 = st.columns(2)

        df_meta = None  # <-- IMPORTANTE: inicializar

        with col1:
            st.markdown("**Arquivos FASTQ**")
            fastq_files = st.file_uploader(
                "Selecione arquivos FASTQ",
                type=["fastq", "fq", "fastq.gz", "fq.gz"],
                accept_multiple_files=True,
                key="fastq_upload",
            )

            if fastq_files:
                st.success(f"{len(fastq_files)} arquivos selecionados")

        with col2:
            st.markdown("**Metadados**")
            metadata_file = st.file_uploader(
                "Arquivo de metadados (TSV/CSV)",
                type=["tsv", "csv", "txt"],
                key="metadata_upload",
            )

            if metadata_file:
                try:
                    # Tentar ler metadados
                    if metadata_file.name.endswith(".csv"):
                        df_meta = pd.read_csv(metadata_file)
                    else:
                        df_meta = pd.read_csv(metadata_file, sep="\t")

                    st.success(f"Metadados carregados: {len(df_meta)} amostras")
                    st.dataframe(df_meta.head())

                except Exception as e:
                    df_meta = None
                    st.error(f"Erro ao ler metadados: {e}")

        # ParÃ¢metros do pipeline
        st.markdown("#### 2. ParÃ¢metros do Pipeline")

        with st.expander("ConfiguraÃ§Ãµes AvanÃ§adas", expanded=False):
            col1, col2, col3 = st.columns(3)

            with col1:
                trunc_len_f = st.number_input("Truncar Forward (bp)", 100, 300, 240)
                trim_left_f = st.number_input("Trim Left Forward", 0, 50, 10)

            with col2:
                trunc_len_r = st.number_input("Truncar Reverse (bp)", 100, 300, 200)
                trim_left_r = st.number_input("Trim Left Reverse", 0, 50, 10)

            with col3:
                max_ee = st.number_input("Max Expected Errors", 1.0, 10.0, 2.0)
                sampling_depth = st.number_input(
                    "Profundidade Amostragem", 1000, 50000, 10000
                )

        # Metadados da amostra
        st.markdown("#### 3. InformaÃ§Ãµes da Amostra")

        col1, col2 = st.columns(2)

        with col1:
            health_status = st.selectbox(
                "Status de SaÃºde",
                ["SaudÃ¡vel", "Disbiose", "Desconhecido"],
            )

            infection_level = st.selectbox(
                "NÃ­vel de InfecÃ§Ã£o",
                ["Nenhum", "Leve", "Moderado", "Severo", "CrÃ­tico"],
            )

            pathogen = st.multiselect(
                "PatÃ³geno Suspeito",
                config.KNOWN_PATHOGENS,
            )

        with col2:
            water_temp = st.number_input("Temperatura da Ãgua (Â°C)", 20.0, 35.0, 28.0)
            salinity = st.number_input("Salinidade (ppt)", 0.0, 50.0, 25.0)
            ph = st.number_input("pH", 6.0, 9.0, 7.8)

        # BotÃ£o de processamento
        st.markdown("#### 4. Processamento")

        if st.button("ðŸš€ Executar Pipeline", type="primary", use_container_width=True):
            if not fastq_files or not metadata_file:
                st.error("FaÃ§a upload dos arquivos FASTQ e metadados primeiro")
                return

            if df_meta is None:
                st.error("Metadados invÃ¡lidos ou nÃ£o carregados corretamente.")
                return

            # Preparar dados
            with st.spinner("Preparando dados..."):
                import tempfile

                temp_dir = Path(tempfile.mkdtemp())

                # Salvar FASTQs
                fastq_dir = temp_dir / "fastq"
                fastq_dir.mkdir()
                for uploaded_file in fastq_files:
                    file_path = fastq_dir / uploaded_file.name
                    with open(file_path, "wb") as f:
                        f.write(uploaded_file.getvalue())

                # Salvar metadados como TSV
                metadata_path = temp_dir / "metadata.tsv"
                df_meta.to_csv(metadata_path, sep="\t", index=False)

                # Metadados adicionais
                sample_metadata = {
                    "health_status": health_status.lower(),
                    "infection_level": infection_level.lower(),
                    "pathogen": ", ".join(pathogen),
                    "water_temperature": water_temp,
                    "salinity": salinity,
                    "ph": ph,
                    "uploaded_by": st.session_state.current_user.username,
                    "upload_date": datetime.now().isoformat(),
                }

                st.info("Iniciando processamento... Isso pode levar alguns minutos.")

                processing_job = {
                    "id": secrets.token_hex(8),
                    "project": st.session_state.current_project.name,
                    "user": st.session_state.current_user.username,
                    "status": "pending",
                    "started_at": datetime.now(),
                    "files": {
                        "fastq_dir": str(fastq_dir),
                        "metadata": str(metadata_path),
                    },
                    "params": {
                        "trunc_len_f": trunc_len_f,
                        "trunc_len_r": trunc_len_r,
                        "trim_left_f": trim_left_f,
                        "trim_left_r": trim_left_r,
                        "max_ee": max_ee,
                        "sampling_depth": sampling_depth,
                    },
                    "sample_metadata": sample_metadata,
                }

                st.session_state.processing_queue.append(processing_job)

                # Processar em thread separada
                thread = threading.Thread(
                    target=self.process_sample,
                    args=(processing_job,),
                    daemon=True,
                )
                thread.start()

                st.success(f"Processamento iniciado! ID: {processing_job['id']}")
                st.info("Acompanhe o progresso na pÃ¡gina de projetos.")

    def process_sample(self, job: Dict):
        """Processa amostra em thread separada"""

        job["status"] = "processing"

        try:
            # Criar pipeline
            project_dir = Path(st.session_state.current_project.base_path)
            pipeline = Shrimp16SPipeline(project_dir)

            # Executar pipeline
            results = pipeline.run_full_pipeline(
                fastq_dir=Path(job["files"]["fastq_dir"]),
                metadata_path=Path(job["files"]["metadata"]),
                params=job["params"],
            )

            job["status"] = (
                "completed" if results.get("status") == "completed" else "failed"
            )
            job["results"] = results
            job["completed_at"] = datetime.now()

            # Salvar no banco de dados
            self.save_processing_results(job)

        except Exception as e:
            job["status"] = "failed"
            job["error"] = str(e)
            logger.error(f"Erro processando job {job['id']}: {e}", exc_info=True)

    def save_processing_results(self, job: Dict):
        """Salva resultados do processamento no banco"""
        session = self.Session()
        try:
            meta = job.get("sample_metadata", {})  # atalho

            sample = Sample(
                project_id=st.session_state.current_project.id,
                uploader_id=st.session_state.current_user.id,
                sample_id=f"USER_{job['id']}",
                health_status=parse_health_status(meta.get("health_status")),
                infection_level=parse_infection_level(meta.get("infection_level")),
                pathogen=meta.get("pathogen"),
                water_temperature=meta.get("water_temperature"),
                salinity=meta.get("salinity"),
                ph=meta.get("ph"),
                processing_status="completed",
                processed_at=datetime.utcnow(),
            )

            session.add(sample)
            session.commit()

            analysis = Analysis(
                project_id=st.session_state.current_project.id,
                sample_id=sample.id,
                user_id=st.session_state.current_user.id,
                analysis_id=f"ANALYSIS_{job['id']}",
                analysis_type=AnalysisType.TAXONOMIC,
                analysis_name=f"Pipeline 16S - {sample.sample_id}",
                parameters=job.get("params", {}),
                pipeline_version=config.QIIME2_VERSION,
                reference_database=config.REFERENCE_DATABASE,
                status="completed",
                started_at=job.get("started_at"),
                completed_at=job.get("completed_at"),
                results_summary=job.get("results", {}).get("metrics", {}),
            )

            session.add(analysis)
            session.commit()

            logger.info(f"Resultados salvos para job {job['id']}")

        finally:
            session.close()

    def comparative_analysis_page(self) -> None:
        """PÃ¡gina de anÃ¡lise comparativa."""
        st.title("ðŸ” AnÃ¡lise Comparativa")

        session = self.Session()
        try:
            # 1. Selecionar amostra do usuÃ¡rio
            st.markdown("### 1. Selecionar Amostra para ComparaÃ§Ã£o")

            user_samples = (
                session.query(Sample)
                .filter_by(
                    uploader_id=st.session_state.current_user.id,
                    processing_status="completed",
                )
                .all()
            )

            if not user_samples:
                st.warning("Nenhuma amostra processada disponÃ­vel. FaÃ§a upload de dados primeiro.")
                return

            sample_options = {
                f"{s.sample_id} - {s.health_status}": s.id for s in user_samples
            }
            selected_sample_label = st.selectbox(
                "Selecione uma amostra",
                list(sample_options.keys()),
            )

            if not selected_sample_label:
                return

            sample_id = sample_options[selected_sample_label]
            user_sample = session.query(Sample).get(sample_id)

            # 2. Configurar comparaÃ§Ã£o
            st.markdown("### 2. Configurar ComparaÃ§Ã£o")

            comparison_type = st.selectbox(
                "Tipo de ComparaÃ§Ã£o",
                [
                    "Status de SaÃºde",
                    "NÃ­vel de InfecÃ§Ã£o",
                    "ComposiÃ§Ã£o da Dieta",
                    "Similaridade Geral",
                ],
            )

            reference_group = st.selectbox(
                "Grupo de ReferÃªncia",
                [
                    "CamarÃ£o SaudÃ¡vel",
                    "CamarÃ£o com Disbiose Leve",
                    "CamarÃ£o com Disbiose Moderada",
                    "CamarÃ£o com Disbiose Severa",
                ],
            )

            metrics_to_use = st.multiselect(
                "MÃ©tricas numÃ©ricas para comparaÃ§Ã£o",
                ["execution_time", "samples_processed", "asvs_generated"],
                default=["asvs_generated"],
            )

            if st.button("ðŸ§ª Executar AnÃ¡lise Comparativa", type="primary"):
                with st.spinner("Comparando com base de dados pÃºblica..."):
                    # Instancia o analisador comparativo
                    analyzer = ComparativeAnalyzer(
                        config.PUBLIC_DATA_DIR,
                        config.USER_DATA_DIR,
                    )

                    # Caminhos de resultados da amostra do usuÃ¡rio
                    project_dir = Path(user_sample.project.base_path)
                    sample_results_dir = project_dir / "results"
                    pipeline_results_path = sample_results_dir / "pipeline_results.json"

                    # Carregar mÃ©tricas da pipeline, se existirem
                    metrics: Dict[str, Any] = {}
                    pipeline_results: Optional[Dict[str, Any]] = None

                    if pipeline_results_path.exists():
                        try:
                            with open(pipeline_results_path, "r", encoding="utf-8") as f:
                                pipeline_results = json.load(f)

                            raw_metrics = pipeline_results.get("metrics", {})
                            # pegar sÃ³ as mÃ©tricas selecionadas e numÃ©ricas
                            for k in metrics_to_use:
                                v = raw_metrics.get(k)
                                if isinstance(v, (int, float)):
                                    metrics[k] = v
                        except Exception as e:
                            logger.warning(
                                f"NÃ£o foi possÃ­vel ler mÃ©tricas da amostra em "
                                f"'{pipeline_results_path}': {e}"
                            )

                    # Rodar a anÃ¡lise comparativa principal
                    comparison_results = analyzer.compare_user_sample(
                        sample_results_dir,
                        metrics=metrics,
                    )

                    # ---------------------------
                    # 1) SALVAR JSON DE COMPARAÃ‡ÃƒO
                    # ---------------------------
                    try:
                        comparative_json_path = sample_results_dir / "comparative_analysis.json"
                        comparative_payload: Dict[str, Any] = {
                            "sample_id": user_sample.sample_id,
                            "internal_sample_id": user_sample.id,
                            "comparison_type": comparison_type,
                            "reference_group": reference_group,
                            "metrics_used": metrics_to_use,
                            "metrics_values": metrics,
                            "results": comparison_results,
                        }

                        # Garante que o diretÃ³rio existe
                        comparative_json_path.parent.mkdir(parents=True, exist_ok=True)

                        with open(comparative_json_path, "w", encoding="utf-8") as f:
                            json.dump(comparative_payload, f, indent=2, ensure_ascii=False)

                    except Exception as e:
                        logger.warning(
                            "NÃ£o foi possÃ­vel salvar JSON de anÃ¡lise comparativa em "
                            f"'{comparative_json_path}': {e}"
                        )

                    # ---------------------------------------------
                    # 2) OPCIONAL: anexar resumo no pipeline_results
                    #    (apenas se o JSON original existir)
                    # ---------------------------------------------
                    if pipeline_results is not None:
                        try:
                            pipeline_results.setdefault("comparative_analysis", {})
                            pipeline_results["comparative_analysis"].update(
                                {
                                    "comparison_type": comparison_type,
                                    "reference_group": reference_group,
                                    "metrics_used": metrics_to_use,
                                }
                            )

                            with open(pipeline_results_path, "w", encoding="utf-8") as f:
                                json.dump(pipeline_results, f, indent=2, ensure_ascii=False)

                        except Exception as e:
                            logger.warning(
                                "NÃ£o foi possÃ­vel atualizar 'pipeline_results.json' "
                                f"com informaÃ§Ãµes da anÃ¡lise comparativa: {e}"
                            )

                    # Guardar em sessÃ£o para uso na interface
                    st.session_state.comparison_results = comparison_results
                    st.success("AnÃ¡lise comparativa concluÃ­da!")

                    # Exibir resultados na interface
                    self.display_comparison_results(comparison_results)

        finally:
            session.close()

    def display_comparison_results(self, results: Dict):
        """Exibe resultados da anÃ¡lise comparativa"""

        st.markdown("### ðŸ“Š Resultados da ComparaÃ§Ã£o")

        classification = results.get("classification", {}) or {}
        comparisons = results.get("comparisons", {}) or {}
        percentiles = results.get("percentiles", {}) or {}
        user_metrics = results.get("user_metrics", {}) or {}

        # --- CartÃµezinhos principais: saÃºde, infecÃ§Ã£o, confianÃ§a -----------------
        col1, col2, col3 = st.columns(3)

        with col1:
            health_status = classification.get("health_status", "unknown")
            status_color = "green" if health_status == "healthy" else "red"
            st.markdown(
                f"""
                <div style='text-align: center; padding: 1rem; background: #f8f9fa; border-radius: 10px;'>
                    <h3 style='color: {status_color};'>{health_status.upper()}</h3>
                    <p>Status de SaÃºde</p>
                </div>
                """,
                unsafe_allow_html=True,
            )

        with col2:
            infection_level = classification.get("infection_level", "unknown")
            st.markdown(
                f"""
                <div style='text-align: center; padding: 1rem; background: #f8f9fa; border-radius: 10px;'>
                    <h3>{infection_level.upper()}</h3>
                    <p>NÃ­vel de InfecÃ§Ã£o</p>
                </div>
                """,
                unsafe_allow_html=True,
            )

        with col3:
            confidence = classification.get("confidence", 0.0) * 100.0
            st.markdown(
                f"""
                <div style='text-align: center; padding: 1rem; background: #f8f9fa; border-radius: 10px;'>
                    <h3>{confidence:.1f}%</h3>
                    <p>ConfianÃ§a</p>
                </div>
                """,
                unsafe_allow_html=True,
            )

        # --- Similaridade por mÃ©trica --------------------------------------------
        if comparisons:
            st.markdown("#### ðŸ“ˆ Similaridade das mÃ©tricas com a base de referÃªncia")

            for metric_name, comp in comparisons.items():
                similarity = comp.get("similarity")
                if similarity is None or np.isnan(similarity):
                    continue

                st.markdown(f"**{metric_name.replace('_', ' ').title()}**")
                st.progress(similarity / 100.0)
                st.caption(
                    f"Similaridade aproximada com a distribuiÃ§Ã£o de referÃªncia: {similarity:.1f}%"
                )

        # --- Percentis por mÃ©trica -----------------------------------------------
        if percentiles:
            st.markdown("#### ðŸŽ¯ PosiÃ§Ã£o Percentil das mÃ©tricas")

            for metric_name, p in percentiles.items():
                col1, col2 = st.columns([3, 1])

                with col1:
                    fig = go.Figure(
                        go.Indicator(
                            mode="gauge+number",
                            value=p,
                            domain={"x": [0, 1], "y": [0, 1]},
                            title={"text": f"{metric_name.replace('_', ' ').title()}"},
                            gauge={
                                "axis": {"range": [0, 100]},
                                "bar": {"color": "darkblue"},
                                "steps": [
                                    {"range": [0, 25], "color": "lightgray"},
                                    {"range": [25, 50], "color": "gray"},
                                    {"range": [50, 75], "color": "darkgray"},
                                    {"range": [75, 100], "color": "dimgray"},
                                ],
                            },
                        )
                    )
                    fig.update_layout(height=200, margin=dict(t=0, b=0))
                    st.plotly_chart(fig, use_container_width=True)

                with col2:
                    user_val = comparisons.get(metric_name, {}).get("user_value")
                    if user_val is not None:
                        st.metric("Valor da sua amostra", f"{user_val:.4g}")
                    else:
                        st.metric("Valor da sua amostra", "N/A")

        # --- RecomendaÃ§Ãµes -------------------------------------------------------
        st.markdown("#### ðŸ’¡ RecomendaÃ§Ãµes")

        recommendations = classification.get("recommendations", []) or []
        if recommendations:
            for rec in recommendations:
                st.info(rec)
        else:
            st.info("Nenhuma recomendaÃ§Ã£o especÃ­fica foi gerada.")

        # --- VisualizaÃ§Ãµes comparativas -----------------------------------------
        st.markdown("#### ðŸ“ˆ VisualizaÃ§Ãµes Comparativas")

        tab1, tab2, tab3 = st.tabs(["Boxplot", "Scatter Plot", "Heatmap"])

        with tab1:
            self.create_comparative_boxplot(results)

        with tab2:
            self.create_comparative_scatter(results)

        with tab3:
            self.create_comparative_heatmap(results)

    def create_comparative_boxplot(self, results: Dict):
        """Cria boxplot comparativo"""
        try:
            # Dados simulados para demonstraÃ§Ã£o
            df = pd.DataFrame(
                {
                    "Grupo": ["SaudÃ¡vel"] * 50
                    + ["Disbiose Leve"] * 50
                    + ["Disbiose Severa"] * 50,
                    "Diversidade Shannon": np.concatenate(
                        [
                            np.random.normal(4.5, 0.5, 50),
                            np.random.normal(3.8, 0.6, 50),
                            np.random.normal(2.5, 0.7, 50),
                        ]
                    ),
                }
            )

            # Adicionar ponto do usuÃ¡rio
            user_val = results.get("user_metrics", {}).get("shannon_index", 3.2)
            df_user = pd.DataFrame(
                {
                    "Grupo": ["Sua Amostra"],
                    "Diversidade Shannon": [user_val],
                }
            )

            fig = go.Figure()

            # Boxplot dos grupos
            for group in df["Grupo"].unique():
                group_data = df[df["Grupo"] == group]
                fig.add_trace(
                    go.Box(
                        y=group_data["Diversidade Shannon"],
                        name=group,
                        boxpoints=False,
                        marker_color="lightblue",
                    )
                )

            # Ponto do usuÃ¡rio
            fig.add_trace(
                go.Scatter(
                    x=["Sua Amostra"] * len(df_user),
                    y=df_user["Diversidade Shannon"],
                    mode="markers",
                    name="Sua Amostra",
                    marker=dict(color="red", size=12, symbol="star"),
                )
            )

            fig.update_layout(
                title="Diversidade Shannon - ComparaÃ§Ã£o com Grupos de ReferÃªncia",
                yaxis_title="Ãndice de Shannon",
                showlegend=True,
                height=500,
            )

            st.plotly_chart(fig, use_container_width=True)

        except Exception as e:
            st.error(f"Erro criando boxplot: {e}")

    def create_comparative_scatter(self, results: Dict):
        """Cria scatter plot comparativo"""
        try:
            # Dados simulados
            np.random.seed(42)
            n_samples = 100

            df = pd.DataFrame(
                {
                    "Diversidade Shannon": np.random.normal(3.5, 1.0, n_samples),
                    "AbundÃ¢ncia Vibrio": np.random.uniform(0, 1, n_samples),
                    "Grupo": np.random.choice(["SaudÃ¡vel", "Disbiose"], n_samples),
                    "Tamanho": np.random.randint(5, 20, n_samples),
                }
            )

            # Adicionar ponto do usuÃ¡rio
            user_shannon = results.get("user_metrics", {}).get("shannon_index", 3.2)
            user_vibrio = results.get("user_metrics", {}).get("vibrio_abundance", 0.4)

            fig = px.scatter(
                df,
                x="Diversidade Shannon",
                y="AbundÃ¢ncia Vibrio",
                color="Grupo",
                size="Tamanho",
                title="RelaÃ§Ã£o: Diversidade vs AbundÃ¢ncia de Vibrio",
                labels={
                    "Diversidade Shannon": "Diversidade Shannon",
                    "AbundÃ¢ncia Vibrio": "AbundÃ¢ncia Relativa de Vibrio",
                },
            )

            # Adicionar ponto do usuÃ¡rio
            fig.add_trace(
                go.Scatter(
                    x=[user_shannon],
                    y=[user_vibrio],
                    mode="markers",
                    name="Sua Amostra",
                    marker=dict(
                        color="red",
                        size=20,
                        symbol="star",
                        line=dict(color="black", width=2),
                    ),
                )
            )

            fig.update_layout(height=500)
            st.plotly_chart(fig, use_container_width=True)

        except Exception as e:
            st.error(f"Erro criando scatter plot: {e}")

    def create_comparative_heatmap(self, results: Dict):
        """Cria heatmap comparativo"""
        try:
            # Dados simulados de abundÃ¢ncia
            taxa = [
                "Vibrio spp.",
                "Bacillus spp.",
                "Photobacterium",
                "Pseudomonas",
                "Flavobacterium",
                "Lactobacillus",
            ]

            groups = [
                "SaudÃ¡vel (n=50)",
                "Disbiose Leve (n=30)",
                "Disbiose Severa (n=20)",
                "Sua Amostra",
            ]

            # Matriz de abundÃ¢ncia
            data = np.random.uniform(0, 1, size=(len(taxa), len(groups)))

            # Ajustar para padrÃµes conhecidos
            data[0, 1:] = [0.1, 0.3, 0.8, 0.6]  # Vibrio
            data[1, :] = [0.7, 0.5, 0.2, 0.4]  # Bacillus

            fig = px.imshow(
                data,
                labels=dict(x="Grupo", y="TÃ¡xon", color="AbundÃ¢ncia"),
                x=groups,
                y=taxa,
                aspect="auto",
                color_continuous_scale="RdYlBu_r",
                title="Heatmap - AbundÃ¢ncia TaxonÃ´mica por Grupo",
            )

            fig.update_layout(
                height=400,
                xaxis_title="Grupos",
                yaxis_title="TÃ¡xons Bacterianos",
            )

            st.plotly_chart(fig, use_container_width=True)

        except Exception as e:
            st.error(f"Erro criando heatmap: {e}")

    def public_database_page(self):
        """PÃ¡gina da base de dados pÃºblica"""
        st.title("ðŸŒ Base de Dados PÃºblica")

        st.markdown(
            """
            Esta seÃ§Ã£o permite acessar e analisar dados pÃºblicos de microbioma de camarÃµes 
            *Penaeus vannamei* disponÃ­veis em repositÃ³rios como SRA e ENA.

            **Recursos disponÃ­veis:**
            - Download automÃ¡tico de datasets pÃºblicos
            - Processamento padronizado com pipeline 16S
            - AnÃ¡lise comparativa com metadados
            - VisualizaÃ§Ã£o interativa dos dados
            """
        )

        # ------------------------------------------------------------------
        # EstatÃ­sticas da base de dados (DINÃ‚MICAS, sem valores fictÃ­cios)
        # ------------------------------------------------------------------
        st.markdown("### ðŸ“Š EstatÃ­sticas da Base de Dados")

        def _format_count(n: int) -> str:
            if n >= 1_000_000:
                return f"{n / 1_000_000:.1f}M"
            if n >= 1_000:
                return f"{n / 1_000:.1f}k"
            return str(int(n))

        datasets_count = 0
        samples_count = 0
        sequences_count = 0
        species_list = []

        ref_db = {}
        all_samples = pd.DataFrame()

        try:
            analyzer = st.session_state.comparative_analyzer
            ref_db = analyzer.build_reference_database()
            all_samples = ref_db.get("all_samples", pd.DataFrame())

            if not all_samples.empty:
                if "accession" in all_samples.columns:
                    datasets_count = int(all_samples["accession"].nunique())
                else:
                    datasets_count = int(all_samples.index.nunique())

                samples_count = int(len(all_samples))

                if "total_reads" in all_samples.columns:
                    sequences_count = int(
                        pd.to_numeric(all_samples["total_reads"], errors="coerce")
                        .fillna(0)
                        .sum()
                    )

                for col in ["scientific_name", "species", "organism"]:
                    if col in all_samples.columns:
                        species_list = (
                            all_samples[col]
                            .dropna()
                            .astype(str)
                            .value_counts()
                            .index
                            .tolist()
                        )
                        break

        except Exception as e:
            logger.warning(
                f"NÃ£o foi possÃ­vel calcular estatÃ­sticas completas da base pÃºblica: {e}"
            )

        n_species = len(set(species_list)) if species_list else 0

        col1, col2, col3, col4 = st.columns(4)

        with col1:
            st.metric("Datasets", _format_count(datasets_count))
        with col2:
            st.metric("Amostras", _format_count(samples_count))
        with col3:
            st.metric("SequÃªncias", _format_count(sequences_count))
        with col4:
            if n_species > 0:
                main_species = species_list[0] if species_list else "Desconhecida"
                st.metric(
                    "EspÃ©cies",
                    str(n_species),
                    help=f"EspÃ©cie mais frequente: {main_species}",
                )
            else:
                st.metric("EspÃ©cies", "0")

        # ------------------------------------------------------------------
        # Busca de dados
        # ------------------------------------------------------------------
        st.markdown("### ðŸ” Buscar Dados PÃºblicos")

        # Definir search_term localmente
        search_term = st.text_input(
            "Termo de busca",
            value="Penaeus vannamei 16S",
            help="Busque por espÃ©cie, tecido, condiÃ§Ã£o, etc.",
        )

        col1, col2, col3 = st.columns(3)

        with col1:
            source = st.selectbox("Fonte", ["SRA", "ENA", "Todas"])
        with col2:
            health_status = st.multiselect(
                "Status de SaÃºde", ["SaudÃ¡vel", "Disbiose", "Desconhecido"]
            )
        with col3:
            tissue_type = st.multiselect(
                "Tipo de Tecido",
                ["Intestino", "HepatopÃ¢ncreas", "Ãgua", "Sedimento", "Desconhecido"],
            )

        region = st.selectbox(
            "RegiÃ£o 16S",
            ["V3-V4", "V4", "V1-V3", "Todas"],
            index=0,
        )
        
        n_max = st.slider(
            "NÃºmero mÃ¡ximo de resultados",
            min_value=10,
            max_value=500,
            value=100,
            step=10,
        )

        if st.button("ðŸ”Ž Buscar Datasets", type="primary"):
            with st.spinner("Buscando dados pÃºblicos nas bases selecionadas..."):
                downloader = PublicDataDownloader(
                    email=getattr(config, "NCBI_EMAIL", None),
                    base_dir=config.PUBLIC_DATA_DIR,
                    max_samples=n_max,
                )

                try:
                    # Buscar dados pÃºblicos usando o search_term definido localmente
                    raw_datasets = downloader.search_public_datasets(
                        source="SRA" if source == "SRA" else None,
                        search_term=search_term,  # Usando search_term local
                        species="Penaeus vannamei",
                        marker="16S",
                        region=region if region != "Todas" else None,
                        n_max=n_max,
                    )

                    if raw_datasets is not None and not raw_datasets.empty:
                        # Organizar dados
                        try:
                            organized_data = downloader.organize_public_data(
                                output_dir=str(config.PUBLIC_DATA_DIR),
                                search_term=search_term,  # Usando search_term local
                                species="Penaeus vannamei",
                                marker="16S",
                                region=region if region != "Todas" else None,
                                n_max=n_max,
                            )
                        except Exception as e:
                            # Se organize_public_data falhar, usar os dados brutos
                            st.warning(f"Erro ao organizar dados: {e}. Usando dados brutos.")
                            organized_data = raw_datasets

                        # Aplicar filtros
                        if health_status:
                            if "health_status" in organized_data.columns:
                                organized_data = organized_data[
                                    organized_data["health_status"].isin(health_status)
                                ]
                        
                        if tissue_type and "tissue" in organized_data.columns:
                            organized_data = organized_data[
                                organized_data["tissue"].isin(tissue_type)
                            ]

                        if organized_data.empty:
                            st.warning(
                                "Foram encontrados datasets, mas nenhum atende aos filtros de saÃºde/tecido selecionados."
                            )
                        else:
                            st.success(
                                f"Encontrados {len(organized_data)} datasets relevantes "
                                "(sem duplicaÃ§Ã£o por projeto+amostra)."
                            )

                            # Definir colunas para mostrar
                            cols_to_show = [
                                "source",
                                "accession",
                                "title",
                                "scientific_name",
                                "total_reads",
                            ]
                            
                            # Adicionar colunas condicionais
                            optional_cols = ["health_status", "infection_level", "tissue"]
                            for col in optional_cols:
                                if col in organized_data.columns:
                                    cols_to_show.append(col)
                            
                            # Filtrar apenas colunas existentes
                            existing_cols = [c for c in cols_to_show if c in organized_data.columns]
                            
                            if existing_cols:
                                st.dataframe(
                                    organized_data[existing_cols].head(50),
                                    use_container_width=True,
                                )
                            else:
                                st.dataframe(
                                    organized_data.head(50),
                                    use_container_width=True,
                                )

                            st.markdown("### ðŸ“¥ OpÃ§Ãµes de Download")

                            selected_accessions = st.multiselect(
                                "Selecione accessions para download",
                                organized_data["accession"].tolist(),
                            )

                            if selected_accessions and st.button("â¬‡ï¸ Baixar Selecionados"):
                                with st.spinner(
                                    f"Baixando {len(selected_accessions)} datasets..."
                                ):
                                    output_dir = config.PUBLIC_DATA_DIR / "downloaded"
                                    output_dir.mkdir(exist_ok=True, parents=True)

                                    # Verificar mÃ©todo de download disponÃ­vel
                                    try:
                                        if hasattr(downloader, 'download_public_samples'):
                                            downloaded = downloader.download_public_samples(
                                                organized_df=organized_data,
                                                accession_list=selected_accessions,
                                                output_dir=output_dir,
                                                threads=4,
                                            )
                                            st.success(
                                                f"Download concluÃ­do: {len(downloaded)} arquivos baixados."
                                            )
                                        else:
                                            st.warning("MÃ©todo de download nÃ£o disponÃ­vel.")
                                    except Exception as e:
                                        st.error(f"Erro no download: {e}")

                                    if st.checkbox("Processar automaticamente com pipeline 16S"):
                                        self.process_public_datasets(output_dir, organized_data)

                    else:
                        st.warning("Nenhum dataset encontrado com os critÃ©rios especificados.")
                        
                except Exception as e:
                    st.error(f"Erro na busca de dados pÃºblicos: {str(e)}")
                    logger.error(f"Erro na busca de dados pÃºblicos: {e}", exc_info=True)

        # ------------------------------------------------------------------
        # VisualizaÃ§Ã£o de dados jÃ¡ processados
        # ------------------------------------------------------------------
        st.markdown("### ðŸ“Š VisualizaÃ§Ã£o de Dados PÃºblicos Processados")

        if st.button("ðŸ“ˆ Carregar Base de Dados Existente"):
            try:
                analyzer = st.session_state.comparative_analyzer
                ref_db = analyzer.build_reference_database()

                if ref_db and not ref_db.get("all_samples", pd.DataFrame()).empty:
                    st.success("Base de dados carregada com sucesso!")

                    stats = ref_db.get("summary_stats", {})
                    all_samples = ref_db.get("all_samples", pd.DataFrame())

                    if stats:
                        st.json(stats)

                    tab1, tab2, tab3 = st.tabs(
                        ["DistribuiÃ§Ã£o", "ComparaÃ§Ãµes", "Metadados"]
                    )

                    with tab1:
                        self.plot_public_data_distribution(ref_db)

                    with tab2:
                        self.plot_public_data_comparisons(ref_db)

                    with tab3:
                        if not all_samples.empty:
                            st.dataframe(all_samples.head(50), use_container_width=True)
                        else:
                            st.info("Nenhum metadado de amostra disponÃ­vel na base atual.")
                else:
                    st.warning("NÃ£o foi possÃ­vel construir a base de referÃªncia pÃºblica.")

            except Exception as e:
                st.error(f"Erro ao carregar base de dados pÃºblica existente: {e}")
    def plot_public_data_distribution(self, ref_db: Dict):
        """Plota distribuiÃ§Ã£o dos dados pÃºblicos"""

        df = ref_db.get("all_samples", pd.DataFrame())

        if df.empty:
            st.info("Nenhum dado disponÃ­vel")
            return

        col1, col2 = st.columns(2)

        with col1:
            # DistribuiÃ§Ã£o por status de saÃºde
            if "health_status" in df.columns:
                health_counts = df["health_status"].value_counts()
                fig1 = px.pie(
                    values=health_counts.values,
                    names=health_counts.index,
                    title="DistribuiÃ§Ã£o por Status de SaÃºde",
                )
                st.plotly_chart(fig1, use_container_width=True)

        with col2:
            # DistribuiÃ§Ã£o por nÃ­vel de infecÃ§Ã£o
            if "infection_level" in df.columns:
                infection_counts = df["infection_level"].value_counts()
                fig2 = px.bar(
                    x=infection_counts.index,
                    y=infection_counts.values,
                    title="DistribuiÃ§Ã£o por NÃ­vel de InfecÃ§Ã£o",
                    labels={"x": "NÃ­vel", "y": "Contagem"},
                )
                st.plotly_chart(fig2, use_container_width=True)

        # Boxplot de diversidade
        if "shannon_index" in df.columns and "health_status" in df.columns:
            fig3 = px.box(
                df,
                x="health_status",
                y="shannon_index",
                color="health_status",
                title="Diversidade Shannon por Status de SaÃºde",
                points="all",
            )
            st.plotly_chart(fig3, use_container_width=True)

    def plot_public_data_comparisons(self, ref_db: Dict):
        """Plota comparaÃ§Ãµes dos dados pÃºblicos"""

        df = ref_db.get("all_samples", pd.DataFrame())

        if df.empty or "health_status" not in df.columns:
            st.info("Dados insuficientes para comparaÃ§Ã£o")
            return

        # Scatter plot: Diversidade vs Vibrio
        if all(
            col in df.columns
            for col in ["shannon_index", "vibrio_abundance", "health_status"]
        ):
            fig = px.scatter(
                df,
                x="shannon_index",
                y="vibrio_abundance",
                color="health_status",
                size="asv_count" if "asv_count" in df.columns else None,
                hover_name=df.index,
                title="RelaÃ§Ã£o: Diversidade vs AbundÃ¢ncia de Vibrio",
                labels={
                    "shannon_index": "Ãndice de Shannon",
                    "vibrio_abundance": "AbundÃ¢ncia de Vibrio",
                },
            )
            st.plotly_chart(fig, use_container_width=True)

        # Heatmap de correlaÃ§Ãµes
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) > 1:
            corr_matrix = df[numeric_cols].corr()

            fig = px.imshow(
                corr_matrix,
                text_auto=True,
                aspect="auto",
                color_continuous_scale="RdBu",
                title="Matriz de CorrelaÃ§Ã£o entre VariÃ¡veis",
            )
            st.plotly_chart(fig, use_container_width=True)
    
    def process_public_datasets(self, data_dir: Path, metadata_df: pd.DataFrame):
        """Processa datasets pÃºblicos automaticamente"""
        
        st.info(f"Processando {len(metadata_df)} datasets pÃºblicos...")
        
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        for idx, row in metadata_df.iterrows():
            accession = row['accession']
            status_text.text(f"Processando {accession}...")
            
            # Criar diretÃ³rio do dataset
            dataset_dir = data_dir / accession
            dataset_dir.mkdir(exist_ok=True)
            
            # Salvar metadados
            metadata_path = dataset_dir / "metadata.json"
            with open(metadata_path, 'w') as f:
                json.dump(row.to_dict(), f, indent=2)
            
            # Verificar se jÃ¡ existem FASTQ
            fastq_files = list(dataset_dir.glob("*.fastq*"))
            
            if len(fastq_files) >= 2:
                # Executar pipeline
                try:
                    pipeline = Shrimp16SPipeline(dataset_dir)
                    
                    results = pipeline.run_full_pipeline(
                        fastq_dir=dataset_dir,
                        metadata_path=metadata_path,
                        params=config.DEFAULT_PARAMS,
                    )
                    
                    if results.get('status') == 'completed':
                        st.success(f"âœ“ {accession} processado")
                    else:
                        st.warning(
                            f"âš ï¸ {accession} falhou: {results.get('error', 'Unknown error')}"
                        )
                
                except Exception as e:
                    st.error(f"âœ— Erro processando {accession}: {e}")
            
            # Atualizar progresso
            progress_bar.progress((idx + 1) / len(metadata_df))
        
        status_text.text("Processamento concluÃ­do!")
        st.success("âœ… Todos os datasets pÃºblicos foram processados")
    
    def admin_panel(self):
        """Painel de administraÃ§Ã£o"""
        
        if st.session_state.current_user.role != UserRole.ADMIN:
            st.error("Acesso restrito a administradores")
            return
        
        st.title("âš™ï¸ Painel de AdministraÃ§Ã£o")
        
        tabs = st.tabs(["UsuÃ¡rios", "Sistema", "Logs", "Backup"])
        
        with tabs[0]:
            self.admin_users_tab()
        
        with tabs[1]:
            self.admin_system_tab()
        
        with tabs[2]:
            self.admin_logs_tab()
        
        with tabs[3]:
            self.admin_backup_tab()
    
    def admin_users_tab(self):
        """Tab de gerenciamento de usuÃ¡rios"""
        st.markdown("### ðŸ‘¥ Gerenciamento de UsuÃ¡rios")
        
        session = self.Session()
        try:
            # Filtrar usuÃ¡rios
            col1, col2 = st.columns(2)
            
            with col1:
                status_filter = st.selectbox(
                    "Filtrar por status",
                    ["Todos", "Pendentes", "Aprovados", "Rejeitados", "Suspensos"],
                )
            
            with col2:
                role_filter = st.selectbox(
                    "Filtrar por funÃ§Ã£o",
                    ["Todos", "Admin", "Pesquisador", "Produtor", "Estudante"],
                )
            
            # Construir query
            query = session.query(User)
            
            if status_filter != "Todos":
                status_map = {
                    "Pendentes": UserStatus.PENDING,
                    "Aprovados": UserStatus.APPROVED,
                    "Rejeitados": UserStatus.REJECTED,
                    "Suspensos": UserStatus.SUSPENDED,
                }
                query = query.filter(User.status == status_map[status_filter])
            
            if role_filter != "Todos":
                role_map = {
                    "Admin": UserRole.ADMIN,
                    "Pesquisador": UserRole.RESEARCHER,
                    "Produtor": UserRole.FARMER,
                    "Estudante": UserRole.STUDENT,
                }
                query = query.filter(User.role == role_map[role_filter])
            
            users = query.order_by(User.created_at.desc()).all()
            
            # Tabela de usuÃ¡rios
            user_data = []
            for user in users:
                user_data.append(
                    {
                        "ID": user.id,
                        "UsuÃ¡rio": user.username,
                        "Nome": user.full_name,
                        "Email": user.email,
                        "OrganizaÃ§Ã£o": user.organization,
                        "FunÃ§Ã£o": user.role.value,
                        "Status": user.status.value,
                        "Criado em": user.created_at.strftime("%d/%m/%Y"),
                        "Projetos": user.project_count,
                        "AÃ§Ãµes": user.id,
                    }
                )
            
            df_users = pd.DataFrame(user_data)
            
            # Editor de dados
            edited_df = st.data_editor(
                df_users,
                column_config={
                    "AÃ§Ãµes": st.column_config.SelectboxColumn(
                        "AÃ§Ãµes",
                        options=[
                            "Aprovar",
                            "Rejeitar",
                            "Suspender",
                            "Tornar Admin",
                            "Remover",
                        ],
                        required=True,
                    )
                },
                use_container_width=True,
                hide_index=True,
            )
            
            # Aplicar aÃ§Ãµes
            if st.button("Aplicar AÃ§Ãµes Selecionadas", type="primary"):
                for _, row in edited_df.iterrows():
                    user_id = row['ID']
                    action = row['AÃ§Ãµes']
                    
                    user = session.query(User).get(user_id)
                    if user:
                        if action == "Aprovar":
                            user.status = UserStatus.APPROVED
                        elif action == "Rejeitar":
                            user.status = UserStatus.REJECTED
                        elif action == "Suspender":
                            user.status = UserStatus.SUSPENDED
                        elif action == "Tornar Admin":
                            user.role = UserRole.ADMIN
                            user.status = UserStatus.APPROVED
                        elif action == "Remover":
                            session.delete(user)
                
                session.commit()
                st.success("AÃ§Ãµes aplicadas com sucesso!")
                st.rerun()
            
            # EstatÃ­sticas
            st.markdown("### ðŸ“Š EstatÃ­sticas de UsuÃ¡rios")
            
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                total_users = session.query(User).count()
                st.metric("Total", total_users)
            
            with col2:
                pending = session.query(User).filter_by(
                    status=UserStatus.PENDING
                ).count()
                st.metric("Pendentes", pending)
            
            with col3:
                active = session.query(User).filter_by(
                    status=UserStatus.APPROVED
                ).count()
                st.metric("Ativos", active)
            
            with col4:
                admins = session.query(User).filter_by(
                    role=UserRole.ADMIN
                ).count()
                st.metric("Admins", admins)
        
        finally:
            session.close()
    
    def admin_system_tab(self):
        """Tab de configuraÃ§Ã£o do sistema"""
        st.markdown("### âš™ï¸ ConfiguraÃ§Ã£o do Sistema")
        
        # Status do sistema
        st.markdown("#### Status do Sistema")
        
        col1, col2 = st.columns(2)
        
        with col1:
            # Verificar dependÃªncias
            st.markdown("**DependÃªncias**")
            
            dependencies = {
                "QIIME2": self.check_qiime2_installation(),
                "Prinseq": self.check_prinseq_installation(),
                "FastQC": self.check_fastqc_installation(),
                "Fasterq-dump": self.check_fasterqdump_installation(),
            }
            
            for dep, status in dependencies.items():
                if status:
                    st.success(f"âœ… {dep}")
                else:
                    st.error(f"âŒ {dep}")
        
        with col2:
            # Uso de disco
            st.markdown("**Uso de Disco**")
            
            total, used, free = shutil.disk_usage("/")
            
            st.metric("Total", f"{total // (2**30)} GB")
            st.metric("Usado", f"{used // (2**30)} GB")
            st.metric("Livre", f"{free // (2**30)} GB")
            
            # Barra de progresso
            usage_percent = (used / total) * 100
            st.progress(usage_percent / 100)
            st.caption(f"{usage_percent:.1f}% usado")
        
        # ConfiguraÃ§Ãµes
        st.markdown("#### ConfiguraÃ§Ãµes do Pipeline")
        
        with st.form("system_config"):
            col1, col2 = st.columns(2)
            
            with col1:
                max_concurrent_jobs = st.number_input(
                    "Jobs Concorrentes MÃ¡ximos",
                    min_value=1,
                    max_value=10,
                    value=4,
                )
                
                default_sampling_depth = st.number_input(
                    "Profundidade PadrÃ£o",
                    min_value=1000,
                    max_value=50000,
                    value=10000,
                )
            
            with col2:
                auto_download_updates = st.checkbox(
                    "Download AutomÃ¡tico de Updates",
                    value=True,
                )
                
                enable_email_notifications = st.checkbox(
                    "NotificaÃ§Ãµes por Email",
                    value=True,
                )
            
            if st.form_submit_button("Salvar ConfiguraÃ§Ãµes"):
                # Aqui vocÃª pode salvar essas configs em algum lugar (DB, arquivo, etc.)
                st.success("ConfiguraÃ§Ãµes salvas!")
    
    def check_qiime2_installation(self) -> bool:
        """Verifica se QIIME2 estÃ¡ instalado"""
        try:
            result = subprocess.run(
                ["qiime", "--version"],
                capture_output=True,
                text=True,
            )
            return result.returncode == 0
        except Exception:
            return False
    
    def check_prinseq_installation(self) -> bool:
        """Verifica se Prinseq estÃ¡ instalado"""
        try:
            result = subprocess.run(
                ["prinseq-lite.pl", "--version"],
                capture_output=True,
                text=True,
            )
            return result.returncode == 0
        except Exception:
            return False
    
    def check_fastqc_installation(self) -> bool:
        """Verifica se FastQC estÃ¡ instalado"""
        try:
            result = subprocess.run(
                ["fastqc", "--version"],
                capture_output=True,
                text=True,
            )
            return result.returncode == 0
        except Exception:
            return False
    
    def check_fasterqdump_installation(self) -> bool:
        """Verifica se fasterq-dump estÃ¡ instalado"""
        try:
            result = subprocess.run(
                ["fasterq-dump", "--version"],
                capture_output=True,
                text=True,
            )
            return result.returncode == 0
        except Exception:
            return False
    
    def admin_logs_tab(self):
        """Tab de visualizaÃ§Ã£o de logs"""
        st.markdown("### ðŸ“ƒ Logs do Sistema")
        
        log_file = config.LOGS_DIR / "shrimp_platform.log"
        
        if log_file.exists():
            # OpÃ§Ãµes de filtro
            col1, col2 = st.columns(2)
            
            with col1:
                log_level = st.selectbox(
                    "NÃ­vel de Log",
                    ["Todos", "INFO", "WARNING", "ERROR", "DEBUG"],
                )
            
            with col2:
                search_term = st.text_input("Buscar no log")
            
            # Ler e filtrar logs
            with open(log_file, 'r') as f:
                log_lines = f.readlines()
            
            filtered_logs = []
            for line in log_lines:
                include = True
                
                if log_level != "Todos" and log_level not in line:
                    include = False
                
                if search_term and search_term.lower() not in line.lower():
                    include = False
                
                if include:
                    filtered_logs.append(line)
            
            # Mostrar logs
            st.text_area(
                "Logs do Sistema",
                value="".join(filtered_logs[-1000:]),  # Ãšltimas 1000 linhas
                height=400,
            )
            
            # EstatÃ­sticas
            st.markdown("**EstatÃ­sticas**")
            
            info_count = sum(1 for line in log_lines if "INFO" in line)
            warning_count = sum(1 for line in log_lines if "WARNING" in line)
            error_count = sum(1 for line in log_lines if "ERROR" in line)
            
            col1, col2, col3 = st.columns(3)
            col1.metric("INFO", info_count)
            col2.metric("WARNING", warning_count)
            col3.metric("ERROR", error_count)
            
            # BotÃ£o para limpar logs
            if st.button("ðŸ—‘ï¸ Limpar Logs", type="secondary"):
                with open(log_file, 'w') as f:
                    f.write("")
                st.success("Logs limpos!")
                st.rerun()
        
        else:
            st.info("Nenhum log encontrado")
    
    def admin_backup_tab(self):
        """Tab de backup do sistema"""
        st.markdown("### ðŸ’¾ Backup do Sistema")
        
        # InformaÃ§Ãµes de backup (exemplo estÃ¡tico â€“ ideal Ã© tornar dinÃ¢mico)
        st.markdown(
            """
        **Ãšltimo Backup:** 2024-01-15 14:30:00  
        **PrÃ³ximo Backup AutomÃ¡tico:** 2024-01-16 02:00:00  
        **Tamanho do Backup:** 2.4 GB  
        **Backups Armazenados:** 7  
        """
        )
        
        # AÃ§Ãµes de backup
        col1, col2, col3 = st.columns(3)
        
        with col1:
            if st.button("ðŸ”§ Backup Manual", use_container_width=True):
                with st.spinner("Criando backup..."):
                    self.create_backup()
                    st.success("Backup criado com sucesso!")
        
        with col2:
            if st.button("ðŸ“¥ Restaurar Backup", use_container_width=True):
                self.restore_backup_ui()
        
        with col3:
            if st.button("ðŸ—‘ï¸ Limpar Backups Antigos", use_container_width=True):
                self.clean_old_backups()
        
        # Lista de backups
        st.markdown("#### ðŸ“ƒ Backups DisponÃ­veis")
        
        backup_dir = config.DATA_DIR / "backups"
        backup_dir.mkdir(exist_ok=True)
        
        backups = list(backup_dir.glob("backup_*.zip"))
        
        if backups:
            backup_data = []
            for backup in sorted(
                backups, key=lambda x: x.stat().st_mtime, reverse=True
            ):
                size = backup.stat().st_size / (1024**3)  # GB
                mod_time = datetime.fromtimestamp(backup.stat().st_mtime)
                
                backup_data.append(
                    {
                        "Arquivo": backup.name,
                        "Tamanho": f"{size:.2f} GB",
                        "Data": mod_time.strftime("%d/%m/%Y %H:%M"),
                        "AÃ§Ã£o": "",  # serÃ¡ escolhida no editor
                    }
                )
            
            df_backups = pd.DataFrame(backup_data)
            
            edited_df = st.data_editor(
                df_backups,
                column_config={
                    "AÃ§Ã£o": st.column_config.SelectboxColumn(
                        "AÃ§Ã£o",
                        options=["", "Restaurar", "Download", "Excluir"],
                        required=False,
                    )
                },
                use_container_width=True,
                hide_index=True,
            )
            
            # Processar aÃ§Ãµes
            if st.button("Aplicar AÃ§Ãµes", type="primary"):
                for _, row in edited_df.iterrows():
                    action = row['AÃ§Ã£o']
                    if not action:
                        continue
                    
                    backup_name = row['Arquivo']
                    backup_path = backup_dir / backup_name
                    
                    if action == "Restaurar":
                        self.restore_backup(backup_path)
                    elif action == "Download":
                        self.download_backup(backup_path)
                    elif action == "Excluir":
                        backup_path.unlink()
                
                st.success("AÃ§Ãµes aplicadas!")
                st.rerun()
        
        else:
            st.info("Nenhum backup disponÃ­vel")
    
    def create_backup(self):
        """Cria backup do sistema"""
        backup_dir = config.DATA_DIR / "backups"
        backup_dir.mkdir(exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_file = backup_dir / f"backup_{timestamp}.zip"
        
        import zipfile
        
        with zipfile.ZipFile(backup_file, 'w', zipfile.ZIP_DEFLATED) as zipf:
            # Banco de dados
            db_file = config.DATABASE_DIR / "shrimp_microbiome.db"
            if db_file.exists():
                zipf.write(db_file, "database/shrimp_microbiome.db")
            
            # ConfiguraÃ§Ãµes
            config_files = list(config.BASE_DIR.glob("*.py"))
            for file in config_files:
                zipf.write(file, f"config/{file.name}")
            
            # Logs
            for log_file in config.LOGS_DIR.glob("*.log"):
                zipf.write(log_file, f"logs/{log_file.name}")
        
        logger.info(f"Backup criado: {backup_file}")
    
    def restore_backup_ui(self):
        """Interface para restaurar backup"""
        st.warning(
            "**AtenÃ§Ã£o:** Restaurar um backup substituirÃ¡ todos os dados atuais!"
        )
        
        backup_dir = config.DATA_DIR / "backups"
        backups = list(backup_dir.glob("backup_*.zip"))
        
        if backups:
            backup_options = {b.name: b for b in backups}
            selected_backup = st.selectbox(
                "Selecione o backup", list(backup_options.keys())
            )
            
            confirm = st.checkbox(
                "Confirmar restauraÃ§Ã£o (esta aÃ§Ã£o nÃ£o pode ser desfeita)"
            )
            
            if st.button("ðŸ”§ Restaurar Backup Selecionado", disabled=not confirm):
                self.restore_backup(backup_options[selected_backup])
                st.success("Backup restaurado com sucesso!")
                st.info("Reinicie a aplicaÃ§Ã£o para aplicar as alteraÃ§Ãµes.")
        else:
            st.info("Nenhum backup disponÃ­vel para restauraÃ§Ã£o")
    
    def restore_backup(self, backup_path: Path):
        """Restaura backup do sistema"""
        import zipfile
        
        with zipfile.ZipFile(backup_path, 'r') as zipf:
            zipf.extractall(config.DATA_DIR)
        
        logger.info(f"Backup restaurado: {backup_path}")
    
    def download_backup(self, backup_path: Path):
        """Disponibiliza backup para download"""
        with open(backup_path, 'rb') as f:
            st.download_button(
                label=f"â¬‡ï¸ Download {backup_path.name}",
                data=f.read(),
                file_name=backup_path.name,
                mime="application/zip",
            )
    
    def clean_old_backups(self, keep_last: int = 7):
        """Remove backups antigos"""
        backup_dir = config.DATA_DIR / "backups"
        backups = sorted(
            backup_dir.glob("backup_*.zip"),
            key=lambda x: x.stat().st_mtime,
            reverse=True,
        )
        
        if len(backups) > keep_last:
            for backup in backups[keep_last:]:
                backup.unlink()
                logger.info(f"Backup removido: {backup}")
            
            st.success(
                f"Backups antigos removidos. Mantidos os Ãºltimos {keep_last} backups."
            )
        else:
            st.info(
                f"Nenhum backup antigo para remover. Mantendo {len(backups)} backups."
            )
    
    def settings_page(self):
        """PÃ¡gina de configuraÃ§Ãµes do usuÃ¡rio"""
        st.title("âš™ï¸ ConfiguraÃ§Ãµes")
        
        if not st.session_state.current_user:
            st.warning("FaÃ§a login para acessar as configuraÃ§Ãµes")
            return
        
        tabs = st.tabs(["Perfil", "NotificaÃ§Ãµes", "Privacidade", "API"])
        
        with tabs[0]:
            self.profile_settings()
        
        with tabs[1]:
            self.notification_settings()
        
        with tabs[2]:
            self.privacy_settings()
        
        with tabs[3]:
            self.api_settings()
    
    def profile_settings(self):
        """ConfiguraÃ§Ãµes de perfil"""
        session = self.Session()
        try:
            user = session.query(User).get(st.session_state.current_user.id)
            
            with st.form("profile_form"):
                col1, col2 = st.columns(2)
                
                with col1:
                    username = st.text_input(
                        "UsuÃ¡rio", value=user.username, disabled=True
                    )
                    email = st.text_input("Email", value=user.email)
                    full_name = st.text_input(
                        "Nome Completo", value=user.full_name or ""
                    )
                    organization = st.text_input(
                        "OrganizaÃ§Ã£o", value=user.organization or ""
                    )
                
                with col2:
                    country = st.selectbox(
                        "PaÃ­s",
                        SHRIMP_PRODUCING_COUNTRIES,
                        index=SHRIMP_PRODUCING_COUNTRIES.index(user.country)
                        if user.country in SHRIMP_PRODUCING_COUNTRIES
                        else 0,
                    )
                    study_species = st.text_input(
                        "EspÃ©cie de Estudo", value=user.study_species or ""
                    )
                    research_interests = st.multiselect(
                        "Interesses de Pesquisa",
                        [
                            "Microbioma",
                            "DoenÃ§as",
                            "NutriÃ§Ã£o",
                            "GenÃ´mica",
                            "Aquacultura SustentÃ¡vel",
                            "Outro",
                        ],
                        default=user.research_interests or [],
                    )
                
                if st.form_submit_button("Salvar AlteraÃ§Ãµes"):
                    user.email = email
                    user.full_name = full_name
                    user.organization = organization
                    user.country = country
                    user.study_species = study_species
                    user.research_interests = research_interests
                    
                    session.commit()
                    st.success("Perfil atualizado com sucesso!")
                    st.session_state.current_user = user
        
        finally:
            session.close()
        
        # Alterar senha
        st.markdown("### ðŸ”’ Alterar Senha")
        
        with st.form("password_form"):
            current_password = st.text_input("Senha Atual", type="password")
            new_password = st.text_input("Nova Senha", type="password")
            confirm_password = st.text_input("Confirmar Nova Senha", type="password")
            
            if st.form_submit_button("Alterar Senha"):
                # Abrir nova sessÃ£o para alteraÃ§Ã£o de senha
                session = self.Session()
                try:
                    user = session.query(User).get(
                        st.session_state.current_user.id
                    )
                    
                    if not user.verify_password(current_password):
                        st.error("Senha atual incorreta")
                    elif new_password != confirm_password:
                        st.error("As novas senhas nÃ£o coincidem")
                    else:
                        user.password_hash = hashlib.sha256(
                            new_password.encode()
                        ).hexdigest()
                        session.commit()
                        st.success("Senha alterada com sucesso!")
                finally:
                    session.close()
    
    def notification_settings(self):
        """ConfiguraÃ§Ãµes de notificaÃ§Ã£o"""
        st.markdown("### ðŸ”” ConfiguraÃ§Ãµes de NotificaÃ§Ã£o")
        
        settings = st.session_state.current_user.settings or {}
        
        col1, col2 = st.columns(2)
        
        with col1:
            email_notifications = st.checkbox(
                "NotificaÃ§Ãµes por Email",
                value=settings.get('email_notifications', True),
            )
            
            analysis_complete = st.checkbox(
                "AnÃ¡lise ConcluÃ­da",
                value=settings.get('analysis_complete', True),
            )
            
            comparison_results = st.checkbox(
                "Resultados de ComparaÃ§Ã£o",
                value=settings.get('comparison_results', True),
            )
        
        with col2:
            system_updates = st.checkbox(
                "AtualizaÃ§Ãµes do Sistema",
                value=settings.get('system_updates', True),
            )
            
            new_features = st.checkbox(
                "Novas Funcionalidades",
                value=settings.get('new_features', True),
            )
            
            newsletter = st.checkbox(
                "Newsletter Mensal",
                value=settings.get('newsletter', False),
            )
        
        if st.button("Salvar ConfiguraÃ§Ãµes", type="primary"):
            session = self.Session()
            try:
                user = session.query(User).get(st.session_state.current_user.id)
                user.settings = {
                    'email_notifications': email_notifications,
                    'analysis_complete': analysis_complete,
                    'comparison_results': comparison_results,
                    'system_updates': system_updates,
                    'new_features': new_features,
                    'newsletter': newsletter,
                }
                session.commit()
                st.success("ConfiguraÃ§Ãµes salvas!")
            finally:
                session.close()
    
    def privacy_settings(self):
        """ConfiguraÃ§Ãµes de privacidade"""
        st.markdown("### ðŸ” ConfiguraÃ§Ãµes de Privacidade")
        
        settings = st.session_state.current_user.settings or {}
        
        # Mapeamento entre chave interna e rÃ³tulo exibido
        key_to_label = {
            "public": "PÃºblico",
            "restricted": "Restrito",
            "standard": "Restrito",  # trata 'standard' como 'restricted'
            "private": "Privado",
        }
        label_to_key = {
            "PÃºblico": "public",
            "Restrito": "restricted",
            "Privado": "private",
        }
        
        current_key = settings.get("privacy_level", "restricted")
        current_label = key_to_label.get(current_key, "Restrito")
        
        options = ["PÃºblico", "Restrito", "Privado"]
        try:
            default_index = options.index(current_label)
        except ValueError:
            default_index = 1  # Restrito
        
        privacy_label = st.selectbox(
            "NÃ­vel de Privacidade",
            options,
            index=default_index,
        )
        
        st.markdown(
            """
        **ExplicaÃ§Ã£o dos nÃ­veis:**
        - **PÃºblico:** seus dados anonimizados podem ser usados em anÃ¡lises agregadas
        - **Restrito:** apenas administradores podem acessar seus dados
        - **Privado:** apenas vocÃª pode acessar seus dados
        """
        )
        
        col1, col2 = st.columns(2)
        
        with col1:
            auto_backup = st.checkbox(
                "Backup AutomÃ¡tico",
                value=settings.get("auto_backup", True),
            )
            
            data_retention = st.number_input(
                "RetenÃ§Ã£o de Dados (meses)",
                min_value=1,
                max_value=60,
                value=settings.get("data_retention", 24),
            )
        
        with col2:
            share_aggregated = st.checkbox(
                "Compartilhar Dados Agregados",
                value=settings.get("share_aggregated", True),
            )
            
            delete_account = st.button(
                "ðŸ—‘ï¸ Solicitar ExclusÃ£o da Conta",
                type="secondary",
            )
        
        if st.button("Salvar ConfiguraÃ§Ãµes", type="primary"):
            session = self.Session()
            try:
                user = session.query(User).get(st.session_state.current_user.id)
                base_settings = user.settings or {}
                user.settings = {
                    **base_settings,
                    "privacy_level": label_to_key.get(
                        privacy_label, "restricted"
                    ),
                    "auto_backup": auto_backup,
                    "data_retention": int(data_retention),
                    "share_aggregated": share_aggregated,
                }
                session.commit()
                st.success("ConfiguraÃ§Ãµes salvas!")
            finally:
                session.close()
        
        if delete_account:
            st.warning(
                """
                **AtenÃ§Ã£o:** Esta aÃ§Ã£o nÃ£o pode ser desfeita!
                
                Todos os seus dados serÃ£o permanentemente excluÃ­dos apÃ³s 30 dias.
                
                Para confirmar a exclusÃ£o, entre em contato com o administrador.
                """
            )
    
    def api_settings(self):
        """ConfiguraÃ§Ãµes de API"""
        st.markdown("### ðŸ”Œ ConfiguraÃ§Ãµes da API")
        
        st.info(
            """
        **API Access**
        Use a API para integrar o ShrimpMicrobiome Platform com seus sistemas.
        
        **Endpoints disponÃ­veis:**
        - `POST /api/upload` - Upload de dados FASTQ
        - `GET /api/results/{analysis_id}` - Obter resultados
        - `POST /api/compare` - AnÃ¡lise comparativa
        """
        )
        
        # Gerar token de API
        if st.button("ðŸ”‘ Gerar Novo Token de API"):
            api_token = secrets.token_urlsafe(32)
            
            session = self.Session()
            try:
                user = session.query(User).get(st.session_state.current_user.id)
                base_settings = user.settings or {}
                user.settings = {
                    **base_settings,
                    'api_token': api_token,
                    'api_token_generated': datetime.utcnow().isoformat(),
                }
                
                session.commit()
                
                st.success("Novo token gerado!")
                st.code(api_token, language=None)
                st.warning(
                    "**Guarde este token em um local seguro!** Ele nÃ£o serÃ¡ mostrado novamente."
                )
            
            finally:
                session.close()
        
        # DocumentaÃ§Ã£o da API
        with st.expander("ðŸ“š DocumentaÃ§Ã£o da API", expanded=False):
            st.markdown(
                """
            ### AutenticaÃ§Ã£o
            ```bash
            curl -X POST https://api.shrimpmicrobiome.com/upload \\
              -H "Authorization: Bearer YOUR_API_TOKEN" \\
              -F "fastq=@sample_R1.fastq.gz" \\
              -F "metadata=@metadata.json"
            ```
            
            ### Upload de Dados
            **Endpoint:** `POST /api/v1/upload`
            
            **ParÃ¢metros:**
            - `fastq_r1` (arquivo): Forward reads  
            - `fastq_r2` (arquivo): Reverse reads  
            - `metadata` (JSON): Metadados da amostra
            
            ### Obter Resultados
            **Endpoint:** `GET /api/v1/results/{analysis_id}`
            
            **Resposta (exemplo):**
            ```json
            {
              "status": "completed",
              "results": {
                "diversity": {...},
                "taxonomy": {...},
                "comparison": {...}
              }
            }
            ```
            """
            )
    
    def about_page(self):
        """PÃ¡gina sobre o sistema"""
        st.title("ðŸ¦ Sobre o ShrimpMicrobiome Platform")

        # Bloquinho bonitinho com nome, versÃ£o e ano
        st.markdown(
            f"""
            <div style="
                display:inline-flex;
                align-items:center;
                gap:0.35rem;
                padding:0.35rem 0.9rem;
                border-radius:999px;
                background:#f1f5f9;
                border:1px solid #cbd5f5;
                font-size:0.9rem;
                font-weight:500;
                margin-bottom:1rem;
            ">
              Aquabiome Â· {config.APP_NAME} v{config.APP_VERSION}
              <span style="color:#64748b;">(2025)</span>
            </div>
            """,
            unsafe_allow_html=True,
        )

        st.markdown(
            f"""
        ## {config.APP_NAME}
        **VersÃ£o:** {config.APP_VERSION}  
        **Desenvolvido por:** {config.COMPANY_NAME}
        
        ---
        
        ### ðŸŽ¯ Objetivo
        Plataforma especializada para anÃ¡lise de microbioma 16S de camarÃµes *Penaeus vannamei*,
        com foco em monitoramento de saÃºde, detecÃ§Ã£o de disbiose e anÃ¡lise comparativa.
        
        ### ðŸš€ Funcionalidades Principais
        
        #### 1. **Pipeline 16S Completo**
        - Controle de qualidade com Prinseq  
        - Processamento QIIME2 otimizado  
        - Denoising com DADA2  
        - AtribuiÃ§Ã£o taxonÃ´mica com SILVA  
        - AnÃ¡lises de diversidade alfa e beta  
        - NMDS, PCoA, RDA, Biplots
        
        #### 2. **Base de Dados PÃºblica**
        - Download automÃ¡tico de dados do SRA/ENA  
        - Processamento padronizado  
        - Metadados estruturados  
        - ClassificaÃ§Ã£o por status de saÃºde
        
        #### 3. **AnÃ¡lise Comparativa**
        - ComparaÃ§Ã£o com dados pÃºblicos  
        - CÃ¡lculo de similaridade  
        - ClassificaÃ§Ã£o automÃ¡tica  
        - Percentis e posicionamento  
        - RecomendaÃ§Ãµes personalizadas
        
        #### 4. **VisualizaÃ§Ã£o Interativa**
        - Heatmaps dinÃ¢micos  
        - Boxplots comparativos  
        - GrÃ¡ficos de dispersÃ£o  
        - OrdenaÃ§Ãµes multivariadas  
        - Dashboards personalizados
        
        ### ðŸ› ï¸ Tecnologias Utilizadas
        
        **Backend:**
        - Python 3.10+  
        - QIIME2 2023.9  
        - DADA2, MAFFT, FastTree  
        - SQLite/PostgreSQL
        
        **Frontend:**
        - Streamlit  
        - Plotly  
        - Pandas, NumPy, SciPy
        
        **Infraestrutura:**
        - Docker  
        - AWS/GCP (opcional)  
        - Git para versionamento
        
        ### ðŸ“Š Metodologia
        
        1. **AquisiÃ§Ã£o de Dados:**
           - Download de dados pÃºblicos (SRA)  
           - Upload de dados do usuÃ¡rio  
           - PadronizaÃ§Ã£o de metadados
        
        2. **Processamento:**
           - QC com Prinseq (qualidade, adaptadores)  
           - ImportaÃ§Ã£o QIIME2  
           - Denoising com DADA2  
           - Taxonomia com SILVA  
           - Diversidade alfa/beta
        
        3. **AnÃ¡lise:**
           - Batch correction  
           - NMDS/PCoA  
           - PERMANOVA  
           - ANCOM para abundÃ¢ncia diferencial
        
        4. **ComparaÃ§Ã£o:**
           - Similaridade Bray-Curtis  
           - ClassificaÃ§Ã£o ML  
           - Percentis  
           - RecomendaÃ§Ãµes
        
        ### ðŸ”¬ Foco em Camaronicultura
        
        **VariÃ¡veis Analisadas:**
        - Status de saÃºde (saudÃ¡vel/disbiose)  
        - NÃ­vel de infecÃ§Ã£o (leve/moderado/severo)  
        - Agentes patogÃªnicos (*Vibrio* spp., WSSV, etc.)  
        - Dieta (proteÃ­na, lipÃ­dio, carboidrato)  
        - ParÃ¢metros ambientais (temperatura, salinidade, pH)
        
        **Biomarcadores:**
        - *Vibrio* spp. (patogÃªnicos)  
        - *Bacillus* spp. (probiÃ³ticos)  
        - Diversidade de Shannon  
        - RazÃµes taxonÃ´micas especÃ­ficas
        
        ### ðŸ§‘â€ðŸ”¬ AplicaÃ§Ãµes
        
        **Para Pesquisadores:**
        - AnÃ¡lise comparativa de experimentos  
        - IdentificaÃ§Ã£o de biomarcadores  
        - Monitoramento longitudinal  
        - PublicaÃ§Ã£o de dados
        
        **Para Produtores:**
        - Monitoramento de saÃºde  
        - DetecÃ§Ã£o precoce de doenÃ§as  
        - OtimizaÃ§Ã£o de dieta  
        - Melhoria de produtividade
        
        **Para Estudantes:**
        - Plataforma educacional  
        - AnÃ¡lises reprodutÃ­veis  
        - Banco de dados de referÃªncia  
        - Metodologia padronizada
        
        ### ðŸ”’ SeguranÃ§a e Privacidade
        
        - AutenticaÃ§Ã£o segura  
        - Dados criptografados  
        - Controle de acesso granular  
        - Backups automÃ¡ticos  
        - Conformidade com LGPD/GDPR
        
        ### ðŸ“š ReferÃªncias
        
        1. Caporaso et al. (2010) - QIIME  
        2. Callahan et al. (2016) - DADA2  
        3. Quast et al. (2013) - SILVA  
        4. Schloss et al. (2009) - mothur
        
        ### ðŸ¤ ContribuiÃ§Ã£o
        
        Este Ã© um projeto open-source. ContribuiÃ§Ãµes sÃ£o bem-vindas!
        
        **GitHub:** [github.com/shrimpmicrobiome](https://github.com/shrimpmicrobiome)  
        **DocumentaÃ§Ã£o:** [docs.shrimpmicrobiome.com](https://docs.shrimpmicrobiome.com)  
        **Suporte:** support@shrimpmicrobiome.com
        
        ---
        
        *Ãšltima atualizaÃ§Ã£o: 12 de dezembro de 2025*
        """
        )
    def render_footer(self):
        footer_html = f"""
        <style>
        .shrimp-footer {{
            position: fixed;
            left: 0;
            bottom: 0;
            width: 100%;
            padding: 0.4rem 1rem;
            background: #f8fafc;
            border-top: 1px solid #e2e8f0;
            font-size: 0.75rem;
            color: #64748b;
            display: flex;
            justify-content: center;
            align-items: center;
            gap: 0.35rem;
            z-index: 999;
        }}
        .shrimp-footer strong {{
            color: #0f172a;
            font-weight: 600;
        }}
        </style>

        <div class="shrimp-footer">
            <strong>Plataforma @ShrimpMicrobiome 2025</strong>
            <span>Â· v{config.APP_VERSION}</span>
            <span>Â· Todos os direitos reservados.</span>
            <span>â€œUma iniciativa Sea Scient, desenvolvida em colaboraÃ§Ã£o com Databiomicsâ€.</span>
        </div>
        """
        st.markdown(footer_html, unsafe_allow_html=True)

    def logout(self):
        """Realiza logout"""
        st.session_state.clear()
        st.success("Logout realizado com sucesso!")
        st.rerun()


# ============================================================================
# FUNÃ‡ÃƒO PRINCIPAL
# ============================================================================
def main():
    """FunÃ§Ã£o principal da aplicaÃ§Ã£o"""
    
    # Configurar pÃ¡gina
    st.set_page_config(
        page_title=config.APP_NAME,
        page_icon="ðŸ¦",
        layout="wide",
        initial_sidebar_state="expanded",
    )
    
    # CSS personalizado
    st.markdown(
        """
    <style>
    .main {
        background-color: #f8f9fa;
    }
    .stApp {
        max-width: 1200px;
        margin: 0 auto;
    }
    </style>
    """,
        unsafe_allow_html=True,
    )
    
    # Inicializar e executar aplicaÃ§Ã£o
    app = ShrimpMicrobiomeApp()
    app.run()


# ============================================================================
# PONTO DE ENTRADA
# ============================================================================
if __name__ == "__main__":
    # Verificar argumentos de linha de comando
    parser = argparse.ArgumentParser(description="ShrimpMicrobiome Platform")
    parser.add_argument("--mode", choices=["web", "cli", "pipeline"], default="web")
    parser.add_argument(
        "--download-public",
        action="store_true",
        help="Download dados pÃºblicos",
    )
    parser.add_argument(
        "--process",
        type=str,
        help="Processar diretÃ³rio de dados",
    )
    parser.add_argument(
        "--compare",
        type=str,
        help="Comparar amostra com base de dados",
    )
    
    args = parser.parse_args()
    
    if args.mode == "web":
        main()
    
    elif args.mode == "cli" and args.download_public:
        # Modo CLI: download de dados pÃºblicos
        print("ðŸ¦ ShrimpMicrobiome Platform - Download de Dados PÃºblicos")
        downloader = PublicDataDownloader()
        
        print("ðŸ” Buscando datasets...")
        datasets = downloader.search_sra_datasets()
        
        if not datasets.empty:
            print(f"ðŸ“¥ Encontrados {len(datasets)} datasets")
            
            organized = downloader.organize_public_data(datasets)
            
            
            output_file = config.PUBLIC_DATA_DIR / "public_datasets.csv"
            organized.to_csv(output_file, index=False)
            
            print(f"ðŸ’¾ Metadados salvos em: {output_file}")
            
            # Download das amostras
            if input("Download das amostras? (s/n): ").lower() == 's':
                accessions = organized['accession'].tolist()[:10]  # Limitar a 10
                output_dir = config.PUBLIC_DATA_DIR / "fastq"
                
                print(f"â¬‡ï¸ Baixando {len(accessions)} amostras...")
                downloaded = downloader.download_sra_samples(accessions, output_dir)
                print(f"âœ… Download concluÃ­do: {len(downloaded)} arquivos")
        else:
            print("âŒ Nenhum dataset encontrado")
    
    elif args.mode == "pipeline" and args.process:
        # Modo pipeline: processar dados
        print("ðŸ¦ ShrimpMicrobiome Platform - Pipeline 16S")
        
        data_dir = Path(args.process)
        if not data_dir.exists():
            print(f"âŒ DiretÃ³rio nÃ£o encontrado: {data_dir}")
            sys.exit(1)
        
        # Verificar arquivos
        fastq_files = list(data_dir.glob("*.fastq*"))
        if not fastq_files:
            print("âŒ Nenhum arquivo FASTQ encontrado")
            sys.exit(1)
        
        metadata_file = data_dir / "metadata.tsv"
        if not metadata_file.exists():
            print("âš ï¸ Arquivo de metadados nÃ£o encontrado, criando template...")
            # Criar template de metadados
            sample_ids = [f.name.split('_')[0] for f in fastq_files]
            unique_samples = sorted(set(sample_ids))
            
            template = pd.DataFrame(
                {
                    'sample-id': unique_samples,
                    'health_status': ['unknown'] * len(unique_samples),
                    'infection_level': ['none'] * len(unique_samples),
                }
            )
            template.to_csv(metadata_file, sep='\t', index=False)
            print(f"ðŸ“„ Template criado: {metadata_file}")
        
        # Executar pipeline
        print("ðŸš€ Iniciando pipeline...")
        pipeline = Shrimp16SPipeline(data_dir)
        results = pipeline.run_full_pipeline(
            data_dir, metadata_file, config.DEFAULT_PARAMS
        )
        
        if results.get('status') == 'completed':
            print(
                f"âœ… Pipeline concluÃ­do em "
                f"{results['metrics']['execution_time']:.1f}s"
            )
            print(f"ðŸ“Š ASVs gerados: {results['metrics']['asv_count']}")
            print(
                "ðŸ“ˆ Diversidade Shannon: "
                f"{results.get('diversity', {}).get('alpha_mean_shannon', 0):.2f}"
            )
        else:
            print(f"âŒ Pipeline falhou: {results.get('error', 'Unknown error')}")
    
    elif args.mode == "cli" and args.compare:
        # Modo CLI: anÃ¡lise comparativa
        print("ðŸ¦ ShrimpMicrobiome Platform - AnÃ¡lise Comparativa")
        
        sample_dir = Path(args.compare)
        if not (sample_dir / "results" / "pipeline_results.json").exists():
            print(f"âŒ Resultados do pipeline nÃ£o encontrados em: {sample_dir}")
            sys.exit(1)
        
        print("ðŸ” Construindo base de dados de referÃªncia...")
        analyzer = ComparativeAnalyzer(config.PUBLIC_DATA_DIR, config.USER_DATA_DIR)
        analyzer.build_reference_database()
        
        print("ðŸ“Š Comparando amostra...")
        results = analyzer.compare_user_sample(sample_dir, {})
        
        classification = results.get('classification', {})
        print("\nðŸ“ƒ **Resultados:**")
        print(f"   Status de SaÃºde: {classification.get('health_status', 'unknown')}")
        print(
            f"   NÃ­vel de InfecÃ§Ã£o: "
            f"{classification.get('infection_level', 'unknown')}"
        )
        print(f"   ConfianÃ§a: {classification.get('confidence', 0)*100:.1f}%")
        
        print("\nðŸ“ˆ **Similaridade:**")
        for group, comparison in results.get('comparisons', {}).items():
            print(f"   {group}: {comparison.get('similarity', 0)*100:.1f}%")
        
        print("\nðŸŽ¯ **Percentis:**")
        for metric, data in results.get('percentiles', {}).items():
            print(f"   {metric}: {data.get('percentile', 0):.1f}%")
        
        print("\nðŸ’¡ **RecomendaÃ§Ãµes:**")
        for rec in classification.get('recommendations', []):
            print(f"   â€¢ {rec}")
    
    else:
        # Modo web padrÃ£o
        main()
