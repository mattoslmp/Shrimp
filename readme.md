# Shrimp 16S Explorer

Aplicação full-stack em Python (Streamlit) para:
- buscar metadados no NCBI/SRA;
- enriquecer com links de FASTQ no ENA;
- baixar arquivos FASTQ;
- fazer análise de qualidade básica;
- preparar e executar pipeline QIIME2.

## Objetivo
Resolver o problema de "não baixar artigos/dados de SRA ou ENA" com um fluxo robusto usando retry, validação de links e fallback para HTTPS em links FTP da ENA.

---

## Funcionalidades implementadas

### Modo User (cliente)
1. **Busca científica parametrizável**
   - Campo de espécie (default: `Penaeus vannamei`)
   - Campo de marcador (default: `16S`)
   - Limite de registros SRA
2. **Coleta de metadados NCBI/SRA**
   - consulta por organismo + marcador;
   - extração de `BioProject`, `BioSample`, título, run accession.
3. **Enriquecimento com ENA**
   - consulta `filereport` para obter `fastq_ftp`, `md5`, layout, plataforma.
4. **Download de FASTQ**
   - converte links `ftp://` e `ftp.sra.ebi.ac.uk/...` para HTTPS quando necessário;
   - baixa múltiplos arquivos por run (pares R1/R2).
5. **Controle de qualidade básico**
   - leitura de FASTQ/FASTQ.GZ;
   - métricas: número de reads amostradas, tamanho médio, qualidade média, taxa Q30.
6. **Exportação de relatórios**
   - CSV de metadados consolidados;
   - CSV de quality report.

### Modo Admin
1. **Autenticação simples por senha**
   - senha via `ADMIN_PASSWORD` (fallback `admin123`).
2. **Preparação automática do QIIME2**
   - geração de `manifest.csv` no formato aceito pelo QIIME2.
3. **Geração dos comandos de pipeline**
   - import (`qiime tools import`),
   - sumário de demux,
   - denoise com DADA2,
   - sumário da feature table e sequências representativas.
4. **Execução opcional via interface**
   - roda os comandos e exibe log/retorno por etapa.

---

## Interface gráfica moderna
A interface usa:
- layout em gradiente oceânico;
- cartões translúcidos (glassmorphism);
- fluxos separados por modo (User/Admin);
- métricas de monitoramento em tempo real no topo.

---

## Estrutura de diretórios

```bash
.
├── app.py
├── readme.md
└── data/
    ├── downloads/
    ├── reports/
    └── qiime2/
```

---

## Instalação

```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install streamlit pandas requests biopython
```

> Para rodar QIIME2, use uma instalação oficial via conda/mamba conforme documentação do projeto QIIME2.

---

## Variáveis de ambiente

```bash
export NCBI_EMAIL="seu_email@dominio.com"
export NCBI_API_KEY="sua_api_key_ncbi"   # opcional, recomendado
export ADMIN_PASSWORD="sua_senha_admin"
```

---

## Execução

```bash
streamlit run app.py
```

Abra o navegador e use:
- **User**: busca/baixa/QC;
- **Admin**: manifest + comandos + execução QIIME2.

---

## Fluxo completo para *Penaeus vannamei* e outras espécies

1. Defina espécie (`Penaeus vannamei` por padrão; pode trocar por qualquer outra).
2. Defina marcador (`16S` por padrão; pode usar outros genes).
3. Rode busca no SRA.
4. Enriqueça com links do ENA.
5. Baixe FASTQ.
6. Rode QC básico.
7. No modo Admin, gere manifest e execute pipeline QIIME2.

---

## Observações importantes

- "Artigos" não são o mesmo que datasets de sequenciamento. O app é focado em **dados brutos/metadados** (SRA/ENA/NCBI).
- A robustez de download depende da disponibilidade remota dos links da ENA/SRA.
- QIIME2 exige ambiente e dependências próprias, normalmente em conda.

---

## Próximas melhorias sugeridas

1. Dashboard de visualização taxonômica pós-QIIME2.
2. Banco SQLite/PostgreSQL para histórico multiusuário.
3. Fila assíncrona para downloads e processamento em lote.
4. Integração com Kraken2/Bracken/MetaPhlAn para comparação de pipelines.
