import sys
import os
import re
import glob
import pandas as pd
import matplotlib
matplotlib.use('Qt5Agg')  # backend compatível com PyQt6
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QLabel, QLineEdit, QSpinBox, 
                             QCheckBox, QPushButton, QTextEdit, QMessageBox, 
                             QTabWidget, QGroupBox, QRadioButton, QComboBox, 
                             QFormLayout, QFileDialog, QTableWidget, QTableWidgetItem,
                             QHeaderView, QButtonGroup, QScrollArea, QFrame, QSizePolicy)
from PyQt6.QtCore import QProcess, Qt
from PyQt6.QtGui import QFont

class MaraudersApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Marauders GenoMap - Dashboard v3.6 🧬")
        self.resize(1200, 950)
        
        self.total_threads = os.cpu_count() or 4
        try:
            mem = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
            self.total_ram = int(mem / (1024**3))
        except: self.total_ram = 16
        
        self.init_ui()

    def init_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QVBoxLayout(central)
        self.tabs = QTabWidget()
        main_layout.addWidget(self.tabs)

        # --- ABA 1: DOWNLOAD (Azul) ---
        self.tab_dl = QWidget()
        dl_lay = QVBoxLayout(self.tab_dl)
        dl_lay.addWidget(QLabel("🆔 ID do SRA (Accession):"))
        self.sra_in = QLineEdit(); self.sra_in.setMinimumHeight(40); dl_lay.addWidget(self.sra_in)

        m_box = QGroupBox("Configuração de Dados"); m_lay = QVBoxLayout()
        m_lay.addWidget(QLabel("Tipo de amostra:"))
        self.combo_type = QComboBox()
        self.combo_type.clear()
        self.combo_type.addItems([
            "Genômica (WGS / Isolados)",
            "Metagenômica",
            "Transcriptômica (RNA-Seq)",
            "Metatranscriptômica (Amostras)"
        ])
        self.combo_type.currentIndexChanged.connect(self.sync_assembly_ui); m_lay.addWidget(self.combo_type)

        self.lay_w = QWidget(); l_lay = QVBoxLayout(self.lay_w); l_lay.addWidget(QLabel("Layout:")); self.combo_layout = QComboBox()
        self.combo_layout.addItems(["Paired-End (R1 + R2)", "Single-End (R1)"]); l_lay.addWidget(self.combo_layout); m_lay.addWidget(self.lay_w); self.lay_w.hide()
        m_box.setLayout(m_lay); dl_lay.addWidget(m_box); dl_lay.addStretch()
        
        btn_dl_lay = QHBoxLayout()
        self.btn_dl = QPushButton("🛰️  1. Baixar Dados")
        self.btn_dl.setStyleSheet("background-color: #2980b9; color: white; font-weight: bold; height: 60px; border-radius: 10px;")
        self.btn_dl.clicked.connect(self.run_dl)
        self.btn_abort = QPushButton("🛑 Abortar"); self.btn_abort.setEnabled(False)
        self.btn_abort.setStyleSheet("background-color: #c0392b; color: white; font-weight: bold; height: 60px; border-radius: 10px;")
        self.btn_abort.clicked.connect(self.abort_process)
        btn_dl_lay.addWidget(self.btn_dl, 3); btn_dl_lay.addWidget(self.btn_abort, 1); dl_lay.addLayout(btn_dl_lay)
        self.tabs.addTab(self.tab_dl, "1. Download")

        # --- ABA 2: ASSEMBLY (Verde) ---
        self.tab_asm = QWidget(); asm_lay = QVBoxLayout(self.tab_asm)
        self.info_lab = QLabel("💡 Dica: Para genomas isolados, use MEGAHIT ou SPAdes."); self.info_lab.setStyleSheet("background: #e1f5fe; color: #01579b; padding: 12px; border-radius: 5px; font-weight: bold; border: 1px solid #b3e5fc;")
        asm_lay.addWidget(self.info_lab)

        self.group_alg = QGroupBox("Algoritmos Disponíveis"); g_lay = QVBoxLayout()
        
        descs = [
            ("MEGAHIT", "<b>MEGAHIT:</b> Projetado para lidar com grandes volumes de dados de forma extremamente eficiente..."),
            ("SPAdes", "<b>SPADES:</b> Ao contrário do MEGAHIT, ele não prioriza a velocidade, mas sim a precisão absoluta..."),
            ("Trinity", "<b>TRINITY:</b> Focado exclusivamente em Transcriptômica e Metatranscriptômica...")
        ]

        for alg, desc in descs:
            h = QHBoxLayout(); rb = QRadioButton(alg)
            if alg == "MEGAHIT": self.radio_mega = rb; rb.setChecked(True)
            elif alg == "SPAdes": self.radio_spa = rb
            else: self.radio_trinity = rb
            lbl_h = QLabel("❓"); lbl_h.setCursor(Qt.CursorShape.WhatsThisCursor); lbl_h.setStyleSheet("color: #2980b9; font-weight: bold;")
            lbl_h.setToolTip(desc); h.addWidget(rb); h.addStretch(); h.addWidget(lbl_h); g_lay.addLayout(h)

        self.group_alg.setLayout(g_lay); asm_lay.addWidget(self.group_alg)
        
        hw_box = QGroupBox("Hardware"); hw_lay = QFormLayout()
        self.spin_t = QSpinBox(); self.spin_t.setRange(2, self.total_threads); self.spin_t.setValue(max(2, int(self.total_threads*0.7)))
        self.spin_r = QSpinBox(); self.spin_r.setRange(4, self.total_ram); self.spin_r.setValue(max(4, int(self.total_ram*0.7)))
        self.hw_advice_lab = QLabel("📊 Razão: 0 GB/Thread"); hw_lay.addRow("Threads:", self.spin_t); hw_lay.addRow("RAM (GB):", self.spin_r); hw_lay.addRow(self.hw_advice_lab)
        hw_box.setLayout(hw_lay); asm_lay.addWidget(hw_box)
        self.spin_t.valueChanged.connect(self.update_hardware_advice); self.spin_r.valueChanged.connect(self.update_hardware_advice)

        self.c_raw = QCheckBox("🧹 Excluir arquivos brutos ao finalizar"); self.c_trim = QCheckBox("🗑️ Excluir pasta de trimming ao finalizar")
        asm_lay.addWidget(self.c_raw); asm_lay.addWidget(self.c_trim); asm_lay.addStretch()
        self.btn_asm = QPushButton("🚀 Iniciar Pipeline"); self.btn_asm.setStyleSheet("background-color: #27ae60; color: white; font-weight: bold; height: 60px; border-radius: 10px;")
        self.btn_asm.clicked.connect(self.run_asm); asm_lay.addWidget(self.btn_asm)
        self.tabs.addTab(self.tab_asm, "2. Assembly")

        # --- ABA 3: PROTEIN SEARCH (Roxo) ---
        self.tab_ps = QWidget(); ps_lay = QVBoxLayout(self.tab_ps)
        ps_lay.addWidget(QLabel("🧬 Perfil HMM (.hmm):"))
        self.path_hmm = QLineEdit(); btn_h = QPushButton("Selecionar"); btn_h.clicked.connect(lambda: self.select_file(self.path_hmm))
        h2 = QHBoxLayout(); h2.addWidget(self.path_hmm); h2.addWidget(btn_h); ps_lay.addLayout(h2)
        l_pf = QLabel('<a href="https://www.ebi.ac.uk/interpro/entry/pfam/">🔗 Link Pfam para arquivos .hmm</a>')
        l_pf.setOpenExternalLinks(True); ps_lay.addWidget(l_pf); ps_lay.addStretch()
        self.btn_ps = QPushButton("🔍 Iniciar ProtSearch (Automático)"); self.btn_ps.setStyleSheet("background-color: #8e44ad; color: white; font-weight: bold; height: 50px; border-radius: 8px;")
        self.btn_ps.clicked.connect(self.run_ps); ps_lay.addWidget(self.btn_ps)
        self.tabs.addTab(self.tab_ps, "3. Protein Search")

        # ============================================================
        # --- ABA 4: RESULTADOS (Laranja) - ATUALIZADA ---
        # ============================================================
        self.tab_sev = QWidget()
        sev_lay = QVBoxLayout(self.tab_sev)
        sev_lay.setSpacing(8)

        # 4.0 — Get Results no TOPO
        box_get = QGroupBox("Obter Resultados do Pipeline")
        box_get_lay = QVBoxLayout()
        btn_result = QPushButton("📥  Get Results  —  Extrair Sequências dos Hits")
        btn_result.setStyleSheet(
            "background-color: #3498db; color: white; font-weight: bold; "
            "height: 46px; border-radius: 8px; font-size: 13px;"
        )
        btn_result.clicked.connect(self.run_get_res)
        box_get_lay.addWidget(btn_result)
        box_get.setLayout(box_get_lay)
        sev_lay.addWidget(box_get)

        # 4.1 — Tabela de Resultados
        box_tab = QGroupBox("Tabela de Resultados (Top 25 Melhores E-values)")
        box_tab_lay = QVBoxLayout()
        
        h_tbl = QHBoxLayout()
        self.path_tbl = QLineEdit()
        self.path_tbl.setPlaceholderText("Caminho para o arquivo .tbl (deixe em branco para auto-detectar pelo SRA)")
        btn_sel_tbl = QPushButton("Selecionar .tbl")
        btn_sel_tbl.clicked.connect(lambda: self.select_file(self.path_tbl, filtro="Arquivos TBL (*.tbl);;Todos os Arquivos (*)"))
        btn_load_tabela = QPushButton("Carregar Tabela")
        btn_load_tabela.setStyleSheet("background-color: #d35400; color: white; font-weight: bold;")
        btn_load_tabela.clicked.connect(self.tabela_resultados)
        
        h_tbl.addWidget(self.path_tbl)
        h_tbl.addWidget(btn_sel_tbl)
        h_tbl.addWidget(btn_load_tabela)
        
        self.tabela = QTableWidget()
        self.tabela.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.tabela.cellClicked.connect(self.transferir_contig)
        self.tabela.setMaximumHeight(220)

        box_tab_lay.addLayout(h_tbl)
        box_tab_lay.addWidget(self.tabela)
        box_tab.setLayout(box_tab_lay)
        sev_lay.addWidget(box_tab)

        # 4.2 — Severus Snap(e): Extrator de Regiões ATUALIZADO
        box_ext = QGroupBox("✂️  Severus Snap(e) — Extrator de Regiões por Descrição")
        box_ext_lay = QVBoxLayout()
        box_ext_lay.setSpacing(6)

        # Seleção do FASTA de hits
        h_fas = QHBoxLayout()
        self.path_fasta_ext = QLineEdit()
        self.path_fasta_ext.setPlaceholderText("Caminho para o arquivo .fasta de hits (ex: lectin_hits.fasta)...")
        btn_sel_fas = QPushButton("Selecionar .fasta")
        btn_sel_fas.clicked.connect(lambda: self.select_file(self.path_fasta_ext, filtro="Arquivos FASTA (*.fasta *.fa *.faa);;Todos os Arquivos (*)"))
        h_fas.addWidget(self.path_fasta_ext)
        h_fas.addWidget(btn_sel_fas)
        box_ext_lay.addLayout(h_fas)

        # Busca por contig
        h_busca = QHBoxLayout()
        lbl_contig = QLabel("Contig alvo:")
        lbl_contig.setFixedWidth(80)
        self.contig_desejado = QLineEdit()
        self.contig_desejado.setPlaceholderText("Nome da contig (clique na tabela acima ou digite manualmente):")
        h_busca.addWidget(lbl_contig)
        h_busca.addWidget(self.contig_desejado)
        box_ext_lay.addLayout(h_busca)

        # Opção de formato de saída
        h_fmt = QHBoxLayout()
        lbl_fmt = QLabel("Formato de saída:")
        lbl_fmt.setFixedWidth(120)
        self.btn_fmt_aa = QRadioButton("Aminoácidos (.faa)")
        self.btn_fmt_nt = QRadioButton("Nucleotídeos (.fna)")
        self.btn_fmt_aa.setChecked(True)
        fmt_group = QButtonGroup(self)
        fmt_group.addButton(self.btn_fmt_aa)
        fmt_group.addButton(self.btn_fmt_nt)
        
        h_fmt.addWidget(lbl_fmt)
        h_fmt.addWidget(self.btn_fmt_aa)
        h_fmt.addWidget(self.btn_fmt_nt)
        h_fmt.addStretch()
        box_ext_lay.addLayout(h_fmt)

        # Informativo sobre extração de região
        info_ext = QLabel(
            "ℹ️  A região extraída é lida da descrição Prodigal (formato: # start # end # strand). "
            "Modo aminoácidos: salva a proteína do .faa na íntegra. "
            "Modo nucleotídeos: transcreve de volta a proteína para CDS (back-transcription) "
            "usando a sequência proteica — nenhum arquivo de contigs adicional é necessário."
        )
        info_ext.setWordWrap(True)
        info_ext.setStyleSheet("color: #555; font-size: 11px; background: #f9f9f9; padding: 6px; border-radius: 4px;")
        box_ext_lay.addWidget(info_ext)

        # Botão principal de corte
        btn_load_fasta = QPushButton("✂️  Buscar e Extrair Região")
        btn_load_fasta.setStyleSheet(
            "background-color: #e67e22; color: white; font-weight: bold; "
            "height: 42px; border-radius: 8px; font-size: 13px;"
        )
        btn_load_fasta.clicked.connect(self.Severus_snapes)
        box_ext_lay.addWidget(btn_load_fasta)

        box_ext.setLayout(box_ext_lay)
        sev_lay.addWidget(box_ext)
        sev_lay.addStretch()

        self.tabs.addTab(self.tab_sev, "4. Resultados")

        # ============================================================
        # --- ABA 5: RELATÓRIO (Nova) ---
        # ============================================================
        self.tab_rel = QWidget()
        rel_lay = QVBoxLayout(self.tab_rel)
        rel_lay.setSpacing(10)

        # Cabeçalho
        lbl_titulo_rel = QLabel("📋  Relatório de Análise")
        lbl_titulo_rel.setStyleSheet(
            "font-size: 16px; font-weight: bold; color: #2c3e50; "
            "padding: 8px; background: #ecf0f1; border-radius: 6px;"
        )
        rel_lay.addWidget(lbl_titulo_rel)

        # Botão de gerar relatório
        btn_gerar_rel = QPushButton("🔄  Gerar / Atualizar Relatório")
        btn_gerar_rel.setStyleSheet(
            "background-color: #16a085; color: white; font-weight: bold; "
            "height: 46px; border-radius: 8px; font-size: 13px;"
        )
        btn_gerar_rel.clicked.connect(self.gerar_relatorio)
        rel_lay.addWidget(btn_gerar_rel)

        # Área de texto do relatório
        box_txt = QGroupBox("Sumário Estatístico")
        box_txt_lay = QVBoxLayout()
        self.rel_texto = QTextEdit()
        self.rel_texto.setReadOnly(True)
        self.rel_texto.setMinimumHeight(130)
        self.rel_texto.setStyleSheet(
            "font-family: 'Courier New', monospace; font-size: 12px; "
            "background: #fdfefe; color: #2c3e50; border: 1px solid #bdc3c7; "
            "border-radius: 4px; padding: 8px;"
        )
        self.rel_texto.setPlaceholderText(
            "Clique em 'Gerar / Atualizar Relatório' para calcular as estatísticas do projeto atual..."
        )
        box_txt_lay.addWidget(self.rel_texto)
        box_txt.setLayout(box_txt_lay)
        rel_lay.addWidget(box_txt)

        # Área de gráficos — três botões, cada um abre em janela separada
        box_graficos = QGroupBox("📊  Gráficos (abre em nova janela)")
        graficos_lay = QHBoxLayout()
        graficos_lay.setSpacing(10)

        btn_graf_contigs = QPushButton("🧱  Gráfico de Contigs")
        btn_graf_contigs.setStyleSheet(
            "background-color: #27ae60; color: white; font-weight: bold; "
            "height: 70px; border-radius: 8px; font-size: 12px;"
        )
        btn_graf_contigs.clicked.connect(self.abrir_grafico_contigs)

        btn_graf_prots = QPushButton("🧬  Gráfico de Proteínas")
        btn_graf_prots.setStyleSheet(
            "background-color: #8e44ad; color: white; font-weight: bold; "
            "height: 70px; border-radius: 8px; font-size: 12px;"
        )
        btn_graf_prots.clicked.connect(self.abrir_grafico_proteinas)

        graficos_lay.addWidget(btn_graf_contigs)
        graficos_lay.addWidget(btn_graf_prots)
        box_graficos.setLayout(graficos_lay)
        rel_lay.addWidget(box_graficos)
        rel_lay.addStretch()

        self.tabs.addTab(self.tab_rel, "5. Relatório")

        # --- LOG INFERIOR ---
        self.log = QTextEdit(); self.log.setReadOnly(True)
        self.log.document().setMaximumBlockCount(1000) 
        self.log.setStyleSheet("background: black; color: #00FF00; font-family: monospace;")
        self.log.setMaximumHeight(160)
        main_layout.addWidget(self.log)
        self.update_hardware_advice()

    # ============================================================
    # FUNÇÕES DA ABA 4 — RESULTADOS (ATUALIZADA)
    # ============================================================

    def run_get_res(self):
        sid = self.sra_in.text().strip()
        if not sid: 
            QMessageBox.warning(self, "Aviso", "Informe o ID do SRA na Aba 1 antes de obter resultados.")
            return
        
        tbls = glob.glob(f"{sid}/02_HMMER_Results/*.tbl")
        faa = f"{sid}/01_Predicted_Proteins/predicted_proteins.faa"
        
        if not tbls or not os.path.exists(faa):
            QMessageBox.warning(self, "Erro", "Resultados HMMER ou proteínas preditas não encontrados.")
            return
        
        self.execute_proc("bash", [os.path.abspath("scripts/get_Seq_results.sh"), os.path.abspath(tbls[0]), os.path.abspath(faa)], os.path.abspath(sid))

    def tabela_resultados(self):
        self.tabela.clear()
        self.tabela.setRowCount(0)
        
        sid = self.sra_in.text().strip()
        caminho_arquivo = ""

        if sid:
            tbls = glob.glob(f"{sid}/02_HMMER_Results/*.tbl")
            if tbls:
                caminho_arquivo = tbls[0]
                self.path_tbl.setText(caminho_arquivo)
                self.log.append(f">>> [INFO] SRA ID detectado: Carregando tabela de {caminho_arquivo}")
            else:
                self.log.append(f">>> [AVISO] ID {sid} presente, mas a pasta de resultados HMMER está vazia.")
        
        if not caminho_arquivo:
            caminho_arquivo = self.path_tbl.text().strip()
            if not caminho_arquivo:
                QMessageBox.warning(self, "Aviso", "Nenhum SRR ID na Aba 1 e nenhum arquivo selecionado manualmente.")
                return
            self.log.append(f">>> [INFO] Usando caminho manual: {caminho_arquivo}")

        if not os.path.exists(caminho_arquivo):
            QMessageBox.critical(self, "Erro", f"Arquivo não encontrado:\n{caminho_arquivo}")
            return

        nomes_oficiais = [
            'target_name', 'target_accession', 'tlen', 'query_name', 'query_accession', 'qlen',
            'full_E-value', 'full_score', 'full_bias', 'domain_num', 'domain_of',
            'c-Evalue', 'i-Evalue', 'domain_score', 'domain_bias',
            'hmm_from', 'hmm_to', 'ali_from', 'ali_to', 'env_from', 'env_to', 'acc',
            'description'
        ]

        try:
            dados = []
            with open(caminho_arquivo, 'r') as arquivo:
                for linha in arquivo:
                    if not linha.startswith('#'):
                        colunas = linha.strip().split(None, 22)
                        if len(colunas) > 0: dados.append(colunas)

            df = pd.DataFrame(dados, columns=nomes_oficiais)
            df["full_E-value"] = pd.to_numeric(df["full_E-value"])
            ntabela = df[df["full_E-value"] < 0.0001].head(25)

            novatabela = ['target_name', "tlen", "full_E-value", "full_score", "full_bias", "description"]
            tabela_final = ntabela[[c for c in novatabela if c in ntabela.columns]]

            self.tabela.setColumnCount(len(tabela_final.columns))
            self.tabela.setHorizontalHeaderLabels(tabela_final.columns)
            self.tabela.setRowCount(len(tabela_final))

            for linha_idx, linha_dados in enumerate(tabela_final.values):
                for col_idx, valor in enumerate(linha_dados):
                    item = QTableWidgetItem(str(valor))
                    item.setTextAlignment(Qt.AlignmentFlag.AlignCenter) 
                    self.tabela.setItem(linha_idx, col_idx, item)
                    
        except Exception as e:
            QMessageBox.critical(self, "Erro", f"Falha ao processar dados: {str(e)}")

    def transferir_contig(self, row, column):
        item_target = self.tabela.item(row, 0)
        if item_target:
            nome_contig = item_target.text()
            self.contig_desejado.setText(nome_contig)
            self.log.append(f">>> [INFO] Alvo '{nome_contig}' enviado para o extrator Severus Snap(e).")

    def _parse_region_from_description(self, description):
        """
        Tenta extrair start, end e strand da descrição no formato Prodigal:
        # start # end # strand (strand: 1 = forward, -1 = reverse)
        Retorna (start_0based, end_1based, strand) ou None se não encontrado.
        """
        # Formato Prodigal: ID # start # end # strand # ...
        partes = description.split('#')
        if len(partes) >= 4:
            try:
                start = int(partes[1].strip())
                end   = int(partes[2].strip())
                strand_raw = partes[3].strip()
                strand = int(strand_raw) if strand_raw in ('1', '-1') else 1
                return (min(start, end) - 1, max(start, end), strand)
            except ValueError:
                pass
        return None

    def Severus_snapes(self):
        search_term = self.contig_desejado.text().strip()
        if not search_term:
            QMessageBox.warning(self, "Aviso", "Clique em uma linha da tabela ou digite o nome da contig!")
            return

        modo_aa = self.btn_fmt_aa.isChecked()

        # -------- Modo Aminoácidos --------
        if modo_aa:
            sid = self.sra_in.text().strip()
            caminho_fasta = ""

            if sid:
                caminho_auto = os.path.join(sid, "01_Predicted_Proteins", "predicted_proteins.faa")
                if os.path.exists(caminho_auto):
                    caminho_fasta = caminho_auto
                    self.path_fasta_ext.setText(caminho_fasta)
                    self.log.append(f">>> [INFO] Usando proteínas automáticas de {sid}")
                else:
                    self.log.append(f">>> [AVISO] Pasta de proteínas para {sid} não encontrada.")

            if not caminho_fasta:
                caminho_fasta = self.path_fasta_ext.text().strip()
                if not caminho_fasta or not os.path.exists(caminho_fasta):
                    QMessageBox.warning(self, "Aviso", "SRR ID ausente ou pasta inválida. Selecione o .fasta manualmente.")
                    return

            output_file, _ = QFileDialog.getSaveFileName(
                self, "Salvar Sequências (Aminoácidos)", "seqs_extraidas.faa",
                "FASTA Proteínas (*.faa);;FASTA (*.fasta)"
            )
            if not output_file: return

            sequencias_salvas = 0
            try:
                with open(caminho_fasta, "r") as infile, open(output_file, "w") as outfile:
                    for record in SeqIO.parse(infile, "fasta"):
                        if search_term in record.id or search_term in record.description:
                            region = self._parse_region_from_description(record.description)
                            if region:
                                inicio, fim, strand = region
                                # Para .faa já é proteína — cortamos apenas se for possível
                                if len(record.seq) >= (fim - inicio):
                                    record_salvar = record[inicio:fim]
                                else:
                                    record_salvar = record
                                    self.log.append(f">>> [INFO] {record.id}: proteína já em tamanho final, salvando na íntegra.")
                                record_salvar.id = record.id + "_aa_region"
                                record_salvar.description = f"[aa:{inicio+1}-{fim}] " + record.description
                            else:
                                record_salvar = record
                                self.log.append(f">>> [INFO] {record.id}: sem coords na descrição, salvando na íntegra.")
                            SeqIO.write(record_salvar, outfile, "fasta")
                            sequencias_salvas += 1

                msg = f"{sequencias_salvas} sequência(s) de aminoácidos salva(s)!\n{output_file}" if sequencias_salvas > 0 else "Contig não localizada no arquivo."
                (QMessageBox.information if sequencias_salvas > 0 else QMessageBox.warning)(
                    self, "Resultado", msg
                )
            except Exception as e:
                QMessageBox.critical(self, "Erro", f"Erro no processamento: {str(e)}")

        # -------- Modo Nucleotídeos (back-transcription a partir do .faa) --------
        else:
            sid = self.sra_in.text().strip()
            caminho_faa = ""

            if sid:
                caminho_auto = os.path.join(sid, "01_Predicted_Proteins", "predicted_proteins.faa")
                if os.path.exists(caminho_auto):
                    caminho_faa = caminho_auto
                    self.path_fasta_ext.setText(caminho_faa)
                    self.log.append(f">>> [INFO] Usando proteínas automáticas de {sid} para back-transcription.")
                else:
                    self.log.append(f">>> [AVISO] Pasta de proteínas para {sid} não encontrada.")

            if not caminho_faa:
                caminho_faa = self.path_fasta_ext.text().strip()
                if not caminho_faa or not os.path.exists(caminho_faa):
                    QMessageBox.warning(self, "Aviso", "SRR ID ausente ou pasta inválida. Selecione o .faa manualmente.")
                    return

            output_file, _ = QFileDialog.getSaveFileName(
                self, "Salvar CDS (Nucleotídeos back-transcribed)", "seqs_extraidas_cds.fna",
                "FASTA Nucleotídeos (*.fna);;FASTA (*.fasta)"
            )
            if not output_file: return

            # Tabela de back-translation: cada aminoácido → codon mais frequente (tabela 11, procarioto)
            CODON_TABLE = {
                'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT', 'C': 'TGT',
                'Q': 'CAA', 'E': 'GAA', 'G': 'GGT', 'H': 'CAT', 'I': 'ATT',
                'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT', 'P': 'CCT',
                'S': 'TCT', 'T': 'ACT', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT',
                '*': 'TAA', 'X': 'NNN'
            }

            def back_transcribe(aa_seq):
                """Converte sequência proteica → CDS sintético (back-translation)."""
                codons = []
                for aa in str(aa_seq).upper():
                    codons.append(CODON_TABLE.get(aa, 'NNN'))
                return Seq(''.join(codons))

            sequencias_salvas = 0
            try:
                with open(caminho_faa, "r") as infile, open(output_file, "w") as outfile:
                    for record in SeqIO.parse(infile, "fasta"):
                        if search_term in record.id or search_term in record.description:
                            region = self._parse_region_from_description(record.description)

                            # Gera o CDS sintético a partir da sequência proteica
                            cds_seq = back_transcribe(record.seq)

                            # Monta a descrição com as coordenadas genômicas originais
                            if region:
                                inicio, fim, strand = region
                                strand_str = "forward" if strand == 1 else "reverse"
                                desc_nt = f"[genomic_coords:{inicio+1}-{fim}|{strand_str}|back-translated] {record.description}"
                            else:
                                desc_nt = f"[back-translated|coords_not_found] {record.description}"
                                self.log.append(f">>> [AVISO] {record.id}: sem coords na descrição; CDS gerado sem posição genômica.")

                            new_record = SeqRecord(
                                cds_seq,
                                id=record.id + "_cds",
                                description=desc_nt
                            )
                            SeqIO.write(new_record, outfile, "fasta")
                            sequencias_salvas += 1
                            self.log.append(f">>> [NT] {record.id}: {len(record.seq)} aa → {len(cds_seq)} bp CDS gerado.")

                msg = (f"{sequencias_salvas} CDS sintético(s) gerado(s) por back-transcription!\n{output_file}"
                       if sequencias_salvas > 0 else "Contig não localizada no arquivo .faa.")
                (QMessageBox.information if sequencias_salvas > 0 else QMessageBox.warning)(
                    self, "Resultado", msg
                )
            except Exception as e:
                QMessageBox.critical(self, "Erro", f"Erro no processamento NT: {str(e)}")

    # ============================================================
    # FUNÇÕES DA ABA 5 — RELATÓRIO (NOVA)
    # ============================================================

    def _coletar_dados_relatorio(self):
        """
        Varre as pastas do projeto e retorna um dicionário com as estatísticas coletadas.
        """
        sid = self.sra_in.text().strip()
        dados = {
            "sid": sid,
            "contigs_num": 0, "contigs_tamanhos": [],
            "proteinas_num": 0, "proteinas_tamanhos": [],
            "hits_num": 0, "hits_evalues": [],
            "erros": []
        }

        if not sid:
            dados["erros"].append("Nenhum SRA ID informado na Aba 1.")
            return dados

        # --- Contigs (Assembly) ---
        contigs_file = None
        for pattern in [
            f"{sid}/05_Assembly_Results/MEGAHIT_*/final.contigs.fa",
            f"{sid}/05_Assembly_Results/SPADES_*/scaffolds.fasta",
            f"{sid}/05_Assembly_Results/SPADES_*/contigs.fasta",
            f"{sid}/05_Assembly_Results/*.Trinity.fasta",
        ]:
            found = glob.glob(pattern)
            if found: contigs_file = found[0]; break

        if contigs_file and os.path.exists(contigs_file):
            try:
                for rec in SeqIO.parse(contigs_file, "fasta"):
                    dados["contigs_tamanhos"].append(len(rec.seq))
                dados["contigs_num"] = len(dados["contigs_tamanhos"])
            except Exception as e:
                dados["erros"].append(f"Erro lendo contigs: {e}")

        # --- Proteínas Preditas ---
        faa_path = os.path.join(sid, "01_Predicted_Proteins", "predicted_proteins.faa")
        if os.path.exists(faa_path):
            try:
                for rec in SeqIO.parse(faa_path, "fasta"):
                    dados["proteinas_tamanhos"].append(len(rec.seq))
                dados["proteinas_num"] = len(dados["proteinas_tamanhos"])
            except Exception as e:
                dados["erros"].append(f"Erro lendo proteínas: {e}")

        # --- Hits HMMER ---
        tbls = glob.glob(f"{sid}/02_HMMER_Results/*.tbl")
        if tbls:
            try:
                with open(tbls[0]) as fh:
                    for linha in fh:
                        if not linha.startswith('#'):
                            cols = linha.strip().split()
                            if len(cols) >= 7:
                                try:
                                    ev = float(cols[6])
                                    dados["hits_evalues"].append(ev)
                                except ValueError: pass
                dados["hits_num"] = len(dados["hits_evalues"])
            except Exception as e:
                dados["erros"].append(f"Erro lendo hits: {e}")

        return dados

    def gerar_relatorio(self):
        dados = self._coletar_dados_relatorio()
        sid = dados["sid"] or "(não informado)"

        linhas = []
        linhas.append(f"{'='*60}")
        linhas.append(f"  RELATÓRIO DE ANÁLISE — Marauders GenoMap v3.6")
        linhas.append(f"{'='*60}")
        linhas.append(f"  Projeto (SRA ID) : {sid}")
        linhas.append(f"{'─'*60}")

        # Contigs
        linhas.append(f"\n  🧱  CONTIGS MONTADOS")
        linhas.append(f"      Número de contigs : {dados['contigs_num']:>10,}")
        if dados['contigs_tamanhos']:
            sz = dados['contigs_tamanhos']
            linhas.append(f"      Menor contig (bp) : {min(sz):>10,}")
            linhas.append(f"      Maior contig (bp) : {max(sz):>10,}")
            linhas.append(f"      Média (bp)        : {sum(sz)/len(sz):>10,.1f}")
            sorted_sz = sorted(sz, reverse=True)
            cum = 0; n50 = 0; total_bp = sum(sz)
            for s in sorted_sz:
                cum += s
                if cum >= total_bp / 2: n50 = s; break
            linhas.append(f"      N50 (bp)          : {n50:>10,}")

        # Proteínas
        linhas.append(f"\n  🧬  PROTEÍNAS PREDITAS")
        linhas.append(f"      Número de proteínas : {dados['proteinas_num']:>8,}")
        if dados['proteinas_tamanhos']:
            sz = dados['proteinas_tamanhos']
            linhas.append(f"      Menor (aa)          : {min(sz):>8,}")
            linhas.append(f"      Maior (aa)          : {max(sz):>8,}")
            linhas.append(f"      Média  (aa)         : {sum(sz)/len(sz):>8,.1f}")

        # Hits
        linhas.append(f"\n  🎯  HITS HMMER")
        linhas.append(f"      Total de hits       : {dados['hits_num']:>8,}")
        if dados['hits_evalues']:
            sig = [e for e in dados['hits_evalues'] if e < 1e-5]
            linhas.append(f"      Hits E-value < 1e-5 : {len(sig):>8,}")
            linhas.append(f"      Melhor E-value      : {min(dados['hits_evalues']):>8.2e}")

        # Erros
        if dados['erros']:
            linhas.append(f"\n  ⚠️  AVISOS / ARQUIVOS NÃO ENCONTRADOS")
            for err in dados['erros']:
                linhas.append(f"      • {err}")

        linhas.append(f"\n{'='*60}")

        self.rel_texto.setPlainText("\n".join(linhas))
        self.log.append(f">>> [RELATÓRIO] Relatório gerado para {sid}.")

        # Armazena para uso nos gráficos
        self._rel_dados = dados

    def abrir_grafico_contigs(self):
        if not hasattr(self, '_rel_dados'):
            QMessageBox.information(self, "Aviso", "Gere o relatório primeiro.")
            return
        dados = self._rel_dados
        if not dados['contigs_tamanhos']:
            QMessageBox.warning(self, "Aviso", "Nenhum dado de contigs encontrado.")
            return

        sz = sorted(dados['contigs_tamanhos'])
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        fig.suptitle(f"Contigs — {dados['sid']}", fontsize=14, fontweight='bold')
        fig.patch.set_facecolor('#f7fdf9')

        # Histograma de tamanhos
        axes[0].hist(sz, bins=min(60, len(sz)//5+5), color='#27ae60', edgecolor='white', alpha=0.85)
        axes[0].set_title("Distribuição de Tamanhos", fontweight='bold')
        axes[0].set_xlabel("Tamanho (bp)")
        axes[0].set_ylabel("Frequência")
        axes[0].set_facecolor('#eafaf1')
        axes[0].grid(axis='y', linestyle='--', alpha=0.6)

        # Curva N50 acumulada
        total_bp = sum(sz)
        cum = 0; n50_mark = None
        cum_vals = []; bp_vals = sorted(sz, reverse=True)
        for s in bp_vals:
            cum += s; cum_vals.append(cum)
            if n50_mark is None and cum >= total_bp / 2:
                n50_mark = (len(cum_vals), s)
        pct = [c / total_bp * 100 for c in cum_vals]
        axes[1].plot(range(1, len(pct)+1), pct, color='#27ae60', linewidth=2)
        if n50_mark:
            axes[1].axvline(n50_mark[0], color='#e74c3c', linestyle='--', label=f'N50 = {n50_mark[1]:,} bp')
            axes[1].legend()
        axes[1].set_title("Cobertura Acumulada (%)", fontweight='bold')
        axes[1].set_xlabel("Número de Contigs (maiores → menores)")
        axes[1].set_ylabel("% do Genoma Montado")
        axes[1].set_facecolor('#eafaf1')
        axes[1].grid(linestyle='--', alpha=0.5)

        plt.tight_layout()
        plt.show()

    def abrir_grafico_proteinas(self):
        if not hasattr(self, '_rel_dados'):
            QMessageBox.information(self, "Aviso", "Gere o relatório primeiro.")
            return
        dados = self._rel_dados
        if not dados['proteinas_tamanhos']:
            QMessageBox.warning(self, "Aviso", "Nenhum dado de proteínas encontrado.")
            return

        sz = dados['proteinas_tamanhos']
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        fig.suptitle(f"Proteínas Preditas — {dados['sid']}", fontsize=14, fontweight='bold')
        fig.patch.set_facecolor('#fdf5ff')

        # Histograma de tamanhos (aa)
        axes[0].hist(sz, bins=min(60, len(sz)//5+5), color='#8e44ad', edgecolor='white', alpha=0.85)
        axes[0].set_title("Distribuição de Tamanhos (aa)", fontweight='bold')
        axes[0].set_xlabel("Tamanho (aminoácidos)")
        axes[0].set_ylabel("Frequência")
        axes[0].set_facecolor('#f5eef8')
        axes[0].grid(axis='y', linestyle='--', alpha=0.6)

        # Distribuição de hits HMMER por E-value (se disponível)
        if dados['hits_evalues']:
            import numpy as np
            log_ev = [-1 * (abs(e) if e == 0 else (len(str(e).split('e-')[-1]) if 'e-' in str(e) else 1))
                      for e in dados['hits_evalues']]
            # usa log10
            ev_arr = [v for v in dados['hits_evalues'] if v > 0]
            if ev_arr:
                import math
                log_vals = [-math.log10(v) for v in ev_arr]
                axes[1].hist(log_vals, bins=30, color='#e74c3c', edgecolor='white', alpha=0.85)
                axes[1].axvline(5, color='#2c3e50', linestyle='--', label='E-value < 1e-5')
                axes[1].set_title("Distribuição de E-values dos Hits", fontweight='bold')
                axes[1].set_xlabel("-log₁₀(E-value)")
                axes[1].set_ylabel("Número de Hits")
                axes[1].legend()
            else:
                axes[1].text(0.5, 0.5, 'E-values = 0\n(hits perfeitos)', ha='center', va='center', fontsize=12)
                axes[1].axis('off')
        else:
            axes[1].text(0.5, 0.5, 'Sem dados\nde hits HMMER', ha='center', va='center',
                         fontsize=13, color='#8e44ad')
            axes[1].axis('off')
        axes[1].set_facecolor('#f5eef8')

        plt.tight_layout()
        plt.show()

    # ============================================================
    # FUNÇÕES ORIGINAIS DO PIPELINE
    # ============================================================

    def update_hardware_advice(self):
        ratio = self.spin_r.value() / self.spin_t.value()
        status = f"📊 Razão: {ratio:.1f} GB/Thread."
        if ratio < 4.0:
            self.hw_advice_lab.setText(f"{status} ⚠️ Abaixo de 4GB/thread: Risco no Assembly.")
            self.hw_advice_lab.setStyleSheet("color: #e67e22; font-weight: bold;")
        else:
            self.hw_advice_lab.setText(f"{status} ✅ Razão segura (4GB+).")
            self.hw_advice_lab.setStyleSheet("color: #27ae60;")

    def sync_assembly_ui(self, i):
        self.lay_w.setVisible(i >= 2)
        self.radio_mega.setEnabled(True)
        self.radio_spa.setEnabled(i == 0) 
        self.radio_trinity.setEnabled(i >= 2) 
        
        if i == 0:
            self.info_lab.setText("💡 Recomendação: MEGAHIT para rapidez ou SPAdes para isolados puros.")
            self.radio_mega.setChecked(True)
        elif i == 1:
            self.info_lab.setText("💡 Recomendação: MEGAHIT é obrigatório para metagenomas complexos.")
            self.radio_mega.setChecked(True)
        elif i == 2:
            self.info_lab.setText("💡 Recomendação: Trinity é o montador ideal para Transcritos (RNA-Seq).")
            self.radio_trinity.setChecked(True)
        elif i == 3:
            self.info_lab.setText("💡 Recomendação: Trinity reconstrói transcritos, mas MEGAHIT pode ser usado para comunidades.")
            self.radio_trinity.setChecked(True)

    def run_dl(self):
        sid = self.sra_in.text().strip()
        if not sid: return
        lay = "paired" if (self.combo_type.currentIndex() != 2 or self.combo_layout.currentIndex() == 0) else "single"
        self.execute_proc("bash", ["scripts/SRAget.sh", sid, lay])

    def run_asm(self):
        sid = self.sra_in.text().strip()
        if not sid: return
        asm = "trinity" if self.radio_trinity.isChecked() else ("megahit" if self.radio_mega.isChecked() else "spades")
        is_p = (self.combo_layout.currentIndex() == 0 if self.combo_type.currentIndex() == 2 else True)
        r1 = f"{sid}_1.fastq.gz" if is_p else f"{sid}.fastq.gz"
        r2 = f"{sid}_2.fastq.gz" if is_p else "none"
        args = [r1, r2, str(self.spin_t.value()), str(self.spin_r.value()), asm,
                "s" if self.c_raw.isChecked() else "n", "s" if self.c_trim.isChecked() else "n"]
        self.execute_proc("bash", [os.path.abspath("scripts/run_assembly.sh")] + args, os.path.abspath(sid))

    def run_ps(self):
        sid = self.sra_in.text().strip()
        hmm = self.path_hmm.text()
        if not sid or not hmm: 
            QMessageBox.warning(self, "Erro", "Informe o ID do SRA e o arquivo HMM.")
            return
        
        mode = "meta" if self.combo_type.currentIndex() in [1, 3] else "single"
        asm_path = os.path.join(sid, "05_Assembly_Results")
        
        contigs_file = None
        if self.radio_mega.isChecked():
            search = glob.glob(f"{asm_path}/MEGAHIT_*/final.contigs.fa")
            if search: contigs_file = search[0]
        elif self.radio_spa.isChecked():
            search = glob.glob(f"{asm_path}/SPADES_*/scaffolds.fasta")
            if not search: search = glob.glob(f"{asm_path}/SPADES_*/contigs.fasta")
            if search: contigs_file = search[0]
        elif self.radio_trinity.isChecked():
            search = glob.glob(f"{asm_path}/*.Trinity.fasta")
            if search: contigs_file = search[0]

        if not contigs_file or not os.path.exists(contigs_file):
            QMessageBox.warning(self, "Erro", "Arquivo de contigs não encontrado para o montador selecionado.")
            return

        self.log.append(f">>> [INFO] Contigs detectados: {os.path.basename(contigs_file)}")
        self.execute_proc("bash", [os.path.abspath("scripts/ProtSearch.sh"), os.path.abspath(contigs_file),
                                   os.path.abspath(hmm), mode], os.path.abspath(sid))

    def execute_proc(self, cmd, args, wd=None):
        self.set_ui_busy(True); self.proc = QProcess()
        if wd: self.proc.setWorkingDirectory(wd)
        self.proc.readyReadStandardOutput.connect(self.update_log)
        self.proc.readyReadStandardError.connect(self.update_log)
        self.proc.finished.connect(lambda: self.set_ui_busy(False))
        self.proc.start(cmd, args)

    def select_file(self, le, filtro="Todos os Arquivos (*)"):
        f, _ = QFileDialog.getOpenFileName(self, "Selecionar", "", filtro)
        if f: le.setText(f)

    def set_ui_busy(self, b):
        self.btn_dl.setEnabled(not b)
        self.btn_asm.setEnabled(not b)
        self.btn_ps.setEnabled(not b)
        self.btn_abort.setEnabled(b)
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor if b else Qt.CursorShape.ArrowCursor)

    def update_log(self):
        out = self.proc.readAllStandardOutput().data().decode(errors='ignore')
        err = self.proc.readAllStandardError().data().decode(errors='ignore')
        clean = (out + err).replace('\r', '\n')
        
        if clean.strip():
            cursor = self.log.textCursor()
            cursor.movePosition(cursor.MoveOperation.End)
            self.log.setTextCursor(cursor)
            self.log.insertPlainText(clean)
            self.log.ensureCursorVisible()

    def abort_process(self):
        if hasattr(self, 'proc'): self.proc.kill(); self.set_ui_busy(False)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = MaraudersApp()
    w.show()
    sys.exit(app.exec())
