import sys
import os
import re
import glob
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QLabel, QLineEdit, QSpinBox, 
                             QCheckBox, QPushButton, QTextEdit, QMessageBox, 
                             QTabWidget, QGroupBox, QRadioButton, QComboBox, 
                             QFormLayout, QFileDialog, QTableWidget, QTableWidgetItem, QHeaderView)
from PyQt6.QtCore import QProcess, Qt

class MaraudersApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Marauders GenoMap - Dashboard v3.5 🧬")
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
            "Genômica (WGS / Isolados)",         # Index 0
            "Metagenômica",                      # Index 1
            "Transcriptômica (RNA-Seq)",         # Index 2
            "Metatranscriptômica (Amostras)"     # Index 3
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

        # --- ABA 4: GET RESULTS (Ciano) ---
        self.tab_res = QWidget(); res_lay = QVBoxLayout(self.tab_res)
        res_lay.addWidget(QLabel("📥 Esta aba coletará automaticamente os resultados HMMER (.tbl) e Proteínas (.faa) dentro da pasta do SRA."))
        res_lay.addStretch()
        self.btn_res = QPushButton("📥 Extrair Sequências (Automático)"); self.btn_res.setStyleSheet("background-color: #16a085; color: white; font-weight: bold; height: 50px; border-radius: 8px;")
        self.btn_res.clicked.connect(self.run_get_res); res_lay.addWidget(self.btn_res)
        self.tabs.addTab(self.tab_res, "4. Get Results")

        # --- ABA 5: SEVERUS SNAP(E) (Laranja) ---
        self.tab_sev = QWidget()
        sev_lay = QVBoxLayout(self.tab_sev)

        # 5.1 Tabela de Resultados
        box_tab = QGroupBox("Tabela de Resultados (Top 25 Melhores E-values)")
        box_tab_lay = QVBoxLayout()
        
        h_tbl = QHBoxLayout()
        self.path_tbl = QLineEdit()
        self.path_tbl.setPlaceholderText("Caminho para o arquivo .tbl (Deixe em branco para auto-detectar pelo SRA)")
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

        box_tab_lay.addLayout(h_tbl)
        box_tab_lay.addWidget(self.tabela)
        box_tab.setLayout(box_tab_lay)
        sev_lay.addWidget(box_tab)

        # 5.2 Extrator Severus Snap(e)
        box_ext = QGroupBox("Severus Snap(e) Sequences - Extrator de Regiões")
        box_ext_lay = QVBoxLayout()
        
        h_fas = QHBoxLayout()
        self.path_fasta_ext = QLineEdit()
        self.path_fasta_ext.setPlaceholderText("Caminho para o arquivo .fasta de hits (ex: lectin_hits.fasta)...")
        btn_sel_fas = QPushButton("Selecionar .fasta")
        btn_sel_fas.clicked.connect(lambda: self.select_file(self.path_fasta_ext, filtro="Arquivos FASTA (*.fasta *.fa *.faa);;Todos os Arquivos (*)"))
        h_fas.addWidget(self.path_fasta_ext)
        h_fas.addWidget(btn_sel_fas)
        
        h_busca = QHBoxLayout()
        self.contig_desejado = QLineEdit()
        self.contig_desejado.setPlaceholderText("Digite o nome da contig desejada para cortar:")
        btn_load_fasta = QPushButton("✂️ Buscar e Cortar")
        btn_load_fasta.setStyleSheet("background-color: #e67e22; color: white; font-weight: bold; padding: 5px;")
        btn_load_fasta.clicked.connect(self.Severus_snapes)
        h_busca.addWidget(self.contig_desejado)
        h_busca.addWidget(btn_load_fasta)
        
        btn_grafico = QPushButton("📊 Gráfico de Tamanhos")
        btn_grafico.setStyleSheet("background-color: #3498db; color: white; font-weight: bold; padding: 5px;")
        btn_grafico.clicked.connect(self.gerar_grafico_tamanhos)
        h_busca.addWidget(btn_grafico)
        
        box_ext_lay.addLayout(h_fas)
        box_ext_lay.addLayout(h_busca)
        box_ext.setLayout(box_ext_lay)
        sev_lay.addWidget(box_ext)

        self.tabs.addTab(self.tab_sev, "5. Análise & Cortes")

        # --- LOG INFERIOR ---
        self.log = QTextEdit(); self.log.setReadOnly(True);
        self.log.document().setMaximumBlockCount(1000) 
        self.log.setStyleSheet("background: black; color: #00FF00; font-family: monospace;"); main_layout.addWidget(self.log)
        self.update_hardware_advice()

    # ==========================================
    # FUNÇÕES DA NOVA ABA (Severus Snap(e))
    # ==========================================
    # ==========================================
    # FUNÇÕES DA ABA 5 (CONECTADAS À ABA 1)
    # ==========================================
    def tabela_resultados(self):
        self.tabela.clear()
        self.tabela.setRowCount(0)
        
        # 1. Tenta pegar o ID da primeira aba
        sid = self.sra_in.text().strip()
        caminho_arquivo = ""

        if sid:
            # Se houver ID, a prioridade é a busca automática na pasta do projeto
            tbls = glob.glob(f"{sid}/02_HMMER_Results/*.tbl")
            if tbls:
                caminho_arquivo = tbls[0]
                self.path_tbl.setText(caminho_arquivo) # Atualiza o campo visualmente
                self.log.append(f">>> [INFO] SRA ID detectado: Carregando tabela de {caminho_arquivo}")
            else:
                self.log.append(f">>> [AVISO] ID {sid} presente, mas a pasta de resultados HMMER está vazia.")
        
        # 2. Caso não tenha ID ou a busca automática falhou, usa o caminho manual
        if not caminho_arquivo:
            caminho_arquivo = self.path_tbl.text().strip()
            if not caminho_arquivo:
                QMessageBox.warning(self, "Aviso", "Nenhum SRR ID na Aba 1 e nenhum arquivo selecionado manualmente.")
                return
            self.log.append(f">>> [INFO] Usando caminho manual: {caminho_arquivo}")

        if not os.path.exists(caminho_arquivo):
            QMessageBox.critical(self, "Erro", f"Arquivo não encontrado:\n{caminho_arquivo}")
            return

        # --- Lógica de leitura do DataFrame (Mantida) ---
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

            # Colunas de exibição
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

    def Severus_snapes(self):
        # 1. Termo de busca
        search_term = self.contig_desejado.text().strip()
        if not search_term:
            QMessageBox.warning(self, "Aviso", "Clique em uma linha da tabela ou digite o nome da contig!")
            return

        # 2. Conexão com o ID da primeira aba para o FASTA
        sid = self.sra_in.text().strip()
        caminho_fasta = ""

        if sid:
            # Tenta localizar as proteínas preditas automaticamente
            caminho_auto = os.path.join(sid, "01_Predicted_Proteins", "predicted_proteins.faa")
            if os.path.exists(caminho_auto):
                caminho_fasta = caminho_auto
                self.path_fasta_ext.setText(caminho_fasta)
                self.log.append(f">>> [INFO] Usando proteínas automáticas de {sid}")
            else:
                self.log.append(f">>> [AVISO] Pasta de proteínas para {sid} não encontrada.")

        # 3. Fallback manual
        if not caminho_fasta:
            caminho_fasta = self.path_fasta_ext.text().strip()
            if not caminho_fasta or not os.path.exists(caminho_fasta):
                QMessageBox.warning(self, "Aviso", "SRR ID ausente ou pasta inválida. Selecione o .fasta manualmente.")
                return

        # 4. Processo de salvamento e corte
        output_file, _ = QFileDialog.getSaveFileName(self, "Salvar Sequências", "seqs_extraidas.fasta", "FASTA (*.fasta)")
        if not output_file: return

        sequencias_salvas = 0 
        try:
            with open(caminho_fasta, "r") as infile, open(output_file, "w") as outfile:
                for record in SeqIO.parse(infile, "fasta"):
                    
                    if search_term in record.id or search_term in record.description:
                        
                        # Por padrão, assumimos que vamos salvar a sequência inteira
                        record_para_salvar = record 
                        
                        partes = record.description.split('#')
                        if len(partes) >= 3:
                            try:
                                c1, c2 = int(partes[1].strip()), int(partes[2].strip())
                                inicio, fim = min(c1, c2) - 1, max(c1, c2)
                                
                                # A TRAVA DE SEGURANÇA: Só corta se a sequência for grande o suficiente!
                                if len(record.seq) >= fim:
                                    # Modo oficial do Biopython de cortar (preserva ID e features corretamente)
                                    record_para_salvar = record[inicio:fim]
                                    record_para_salvar.id = record.id + "_cortada"
                                    record_para_salvar.description = f"[{inicio+1}-{fim}] " + record.description
                                else:
                                    # Se for menor, já é a proteína traduzida/cortada
                                    self.log.append(f">>> [INFO] {record.id} já é a proteína final. Salvando na íntegra.")
                                    
                            except ValueError: 
                                pass
                        
                        # Escreve o resultado (cortado ou inteiro) no arquivo final
                        SeqIO.write(record_para_salvar, outfile, "fasta")
                        sequencias_salvas += 1 
            
            if sequencias_salvas > 0:
                QMessageBox.information(self, "Sucesso", f"{sequencias_salvas} sequência(s) processada(s) e salva(s)!\n{output_file}")
            else:
                QMessageBox.warning(self, "Não encontrado", "A contig não foi localizada no arquivo selecionado.")
        except Exception as e:
            QMessageBox.critical(self, "Erro", f"Erro no processamento: {str(e)}")
    
    
    def transferir_contig(self, row, column):
        """Pega o nome da contig da linha clicada e envia para o Severus Snap"""
        # A coluna 0 é onde fica o 'target_name'. 
        # Independentemente de qual coluna o usuário clicar, pegamos o item da coluna 0 daquela linha.
        item_target = self.tabela.item(row, 0)
        
        if item_target:
            nome_contig = item_target.text()
            
            # Preenche a caixa de texto lá embaixo
            self.contig_desejado.setText(nome_contig)
            
            # Opcional: Avisa no log preto que a transferência foi feita
            self.log.append(f">>> [INFO] Alvo '{nome_contig}' enviado para o extrator Severus Snap(e).")
    
    def gerar_grafico_tamanhos(self):
        # 1. Pega o caminho do arquivo FASTA da interface
        caminho_fasta = self.path_fasta_ext.text().strip()
        
        # Se estiver vazio, tenta achar automaticamente pelo SRA (mesma lógica do corte)
        if not caminho_fasta:
            sid = self.sra_in.text().strip()
            if sid:
                caminho_auto = os.path.join(sid, "01_Predicted_Proteins", "predicted_proteins.faa")
                if os.path.exists(caminho_auto):
                    caminho_fasta = caminho_auto
                    self.path_fasta_ext.setText(caminho_fasta)

        # Se continuar vazio ou não existir, avisa o usuário
        if not caminho_fasta or not os.path.exists(caminho_fasta):
            QMessageBox.warning(self, "Aviso", "Por favor, selecione um arquivo .fasta válido para gerar o gráfico.")
            return

        # 2. Extrai os tamanhos usando Biopython
        tamanhos = []
        try:
            with open(caminho_fasta, "r") as infile:
                for record in SeqIO.parse(infile, "fasta"):
                    # Aqui é onde a mágica acontece: len(record.seq) pega o tamanho!
                    tamanhos.append(len(record.seq))
            
            if not tamanhos:
                QMessageBox.warning(self, "Aviso", "O arquivo FASTA parece estar vazio.")
                return

            # 3. Gera o gráfico com Matplotlib
            plt.figure(figsize=(10, 6)) # Define o tamanho da janela do gráfico
            
            # Cria um histograma (bins=50 divide os dados em 50 colunas)
            plt.hist(tamanhos, bins=50, color='#9b59b6', edgecolor='black', alpha=0.7)
            
            # Títulos e labels
            plt.title(f"Distribuição de Tamanhos das Sequências\n({os.path.basename(caminho_fasta)})", fontsize=14)
            plt.xlabel("Tamanho da Sequência (pb ou aa)", fontsize=12)
            plt.ylabel("Frequência (Quantidade de Sequências)", fontsize=12)
            plt.grid(axis='y', linestyle='--', alpha=0.7) # Coloca umas linhas de fundo para facilitar a leitura
            
            # Mostra estatísticas básicas no console preto do Marauders
            media = sum(tamanhos) / len(tamanhos)
            maior = max(tamanhos)
            menor = min(tamanhos)
            self.log.append(f">>> [GRÁFICO] Total de seqs: {len(tamanhos)} | Menor: {menor} | Maior: {maior} | Média: {media:.1f}")

            # Abre a janela do gráfico
            plt.show()

        except Exception as e:
            QMessageBox.critical(self, "Erro", f"Ocorreu um erro ao gerar o gráfico:\n{str(e)}")

    
    
    # ==========================================
    # FUNÇÕES ORIGINAIS DO PIPELINE
    # ==========================================
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
        args = [r1, r2, str(self.spin_t.value()), str(self.spin_r.value()), asm, "s" if self.c_raw.isChecked() else "n", "s" if self.c_trim.isChecked() else "n"]
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
        self.execute_proc("bash", [os.path.abspath("scripts/ProtSearch.sh"), os.path.abspath(contigs_file), os.path.abspath(hmm), mode], os.path.abspath(sid))
        
    def run_get_res(self):
        sid = self.sra_in.text().strip()
        if not sid: return
        
        tbls = glob.glob(f"{sid}/02_HMMER_Results/*.tbl")
        faa = f"{sid}/01_Predicted_Proteins/predicted_proteins.faa"
        
        if not tbls or not os.path.exists(faa):
            QMessageBox.warning(self, "Erro", "Resultados HMMER ou proteínas preditas não encontrados."); return
        
        self.execute_proc("bash", [os.path.abspath("scripts/get_Seq_results.sh"), os.path.abspath(tbls[0]), os.path.abspath(faa)], os.path.abspath(sid))

    def execute_proc(self, cmd, args, wd=None):
        self.set_ui_busy(True); self.proc = QProcess()
        if wd: self.proc.setWorkingDirectory(wd)
        self.proc.readyReadStandardOutput.connect(self.update_log); self.proc.readyReadStandardError.connect(self.update_log)
        self.proc.finished.connect(lambda: self.set_ui_busy(False)); self.proc.start(cmd, args)

    def select_file(self, le, filtro="Todos os Arquivos (*)"):
        f, _ = QFileDialog.getOpenFileName(self, "Selecionar", "", filtro)
        if f: le.setText(f)

    def set_ui_busy(self, b):
        self.btn_dl.setEnabled(not b); self.btn_asm.setEnabled(not b); self.btn_ps.setEnabled(not b); self.btn_res.setEnabled(not b); self.btn_abort.setEnabled(b)
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
    app = QApplication(sys.argv); w = MaraudersApp(); w.show(); sys.exit(app.exec())
