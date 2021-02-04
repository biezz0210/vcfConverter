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
        self.textEdit.setText(file_name[0])
        _file.append(file_name[0])

        self.progress_label.setText('selected')
        self.progressBar.setValue(0)    


    def button_convert(self):

        ver = self.comboBox.currentText()

        self.progressBar.setRange(0,0)
        self.progress_label.setText('ongoing')

        try:
            self.progressConvert(_file[-1], ver)
            self.progressBar.setRange(0,100)
            self.progressBar.setValue(100)
            txt = emoji.emojize('Done :thumbs_up:')
            self.progress_label.setText('%s'%txt)

        except:
            self.progress_label.setText('Error: This file is not supported')


    def progressConvert(self, fileName, ver):
        
        df = pd.read_csv('%s'%fileName, compression='gzip', sep='\t',header=None ,error_bad_lines=False,warn_bad_lines=False)
        x,y = df.shape
        df = pd.read_csv('%s'%fileName, compression='gzip', sep='\t',skiprows=x)

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
        def sorted_chr(x):
            try:
                return int(x.split('Gm')[1])
            except:
                return int(x.split('_')[1]) + 999
        df['sort'] = list(map(sorted_chr, df['#CHROM']))
        df_sort = df.sort_values(by=['sort', 'POS'])

        df_sort = df_sort[df_sort.columns[:-1]]

        header = '''##reference=glyma.Wm82.gnm4.4PTR.genome_main.fna\n'''


        name = QFileDialog.getSaveFileName(self, 'Open file', './')
        if name[0] == '':
            name = './sample.vcf'
        else:
            name = name[0] + '.vcf'

        with open(name, 'w') as vcf:
            vcf.write(header)

        df_sort.to_csv(name, sep='\t', mode='a', index=False)




if __name__ == '__main__':
    app = QApplication(sys.argv)
    Dialog = MainDialog()
    Dialog.show()
    app.exec_()
