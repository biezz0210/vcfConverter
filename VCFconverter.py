import sys
from PyQt5.QtWidgets import *
from PyQt5 import uic
from PyQt5.QtCore import pyqtSlot, Qt
import pandas as pd
import emoji

test_ui = r'./VCFconverter.ui'


class MainDialog(QDialog):
    def __init__(self):
        QDialog.__init__(self, None)
        uic.loadUi(test_ui, self)


        global _file
        _file = []

        self.file_open.clicked.connect(self.button_click_open)
        self.convert.clicked.connect(self.button_convert)


    def button_click_open(self):
        file_name = QFileDialog.getOpenFileName(self)
        # exist = self.textEdit.toPlainText()
        self.textEdit.setText(file_name[0])
        # out = 'PATH :  ' + '/'.join(file_name[0].split('/')[:-1]) + '/'
        # self.label_output.setText(out)
        _file.append(file_name[0])


    def button_convert(self):

        ver = self.comboBox.currentText()

        df = pd.read_csv('%s'%_file[-1], compression='gzip', sep='\t',header=None ,error_bad_lines=False,warn_bad_lines=False)
        x,y = df.shape
        df = pd.read_csv('%s'%_file[-1], compression='gzip', sep='\t',skiprows=x)
        
        m = [x.split(';')[0] != 'INDEL' for x in df['INFO']]
        df_snp = df[m]
        
        chrom = [x.split('.')[-1] for x in df_snp['#CHROM']]
        pos = list(map(str,list(df_snp['POS'])))

        def sum_data(x, y):
            return x + '-' + y

        df_snp['ix'] = list(map(sum_data, chrom, pos))

        if ver[-3:] == '1.0':
            df_db = pd.read_csv('./data/Ver1DB.csv')
        else:
            df_db = pd.read_csv('./data/Ver2DB.csv')
        common_ix = set(df_db['POSpre']) & set(df_snp['ix'])
        df_snp_ix = df_snp.set_index('ix')
        df_db_ix = df_db.set_index('POSpre')
        df_db_common = df_db_ix.loc[common_ix]

        chr_fix = list(map(lambda x : x.split('-')[0], df_db_common['POSv4']))
        pos_fix = list(map(lambda x : int(x.split('-')[1]), df_db_common['POSv4']))
        ref = list(df_db_common['REF'])
        alt = list(df_db_common['ALT'])

        df_snp_common = df_snp_ix.loc[common_ix]

        df_snp_common['#CHROM'] = chr_fix
        df_snp_common['POS'] = pos_fix
        df_snp_common['REF'] = ref
        df_snp_common['ALT'] = alt
        
        df = df_snp_common.reset_index()[df_snp_common.reset_index().columns[1:]]
        df_sort = df.sort_values(by=['#CHROM', 'POS'])

        header = '''##reference=glyma.Wm82.gnm4.4PTR.genome_main.fna\n'''

        txt = emoji.emojize('Done :thumbs_up:')
        self.label_result.setText('%s'%txt)


        name = QFileDialog.getSaveFileName(self, 'Open file', './')
        name = name[0] + '.vcf'

        with open(name, 'w') as vcf:
            vcf.write(header)

        df.to_csv(name, sep='\t', mode='a', index=False)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    Dialog = MainDialog()
    Dialog.show()
    app.exec_()
