import sys
from PyQt5.QtWidgets import (QApplication, QLineEdit, QWidget, QToolTip, QPushButton, QMessageBox, QDesktopWidget, 
QMainWindow, qApp, QAction, QMenu, QTextBrowser, QVBoxLayout, QListWidget, QListWidgetItem, QLabel, QPlainTextEdit, QScrollArea, QScrollBar)
from PyQt5.QtGui import QFont, QIcon
from PyQt5.QtCore import pyqtSlot, Qt
import os
import math
import json


rootDir = '../../Desktop/database/'

## a class which will contain all info about species
class Species:
        def __init__(self, data):
                self.taxon = data['taxon'].replace("_", " ")
                self.accession = data['accession']
                self.seqs = []
                for seq in data['sequences']:
                        newSeq = Sequence(seq)
                        self.seqs.append(newSeq)

# add a class for hits from HMM hits
class Sequence:
        def __init__(self, sequenceData):
                self.title = sequenceData['title']
                self.type = sequenceData['type']
                self.pseq = sequenceData['p_seq']
                self.nseq = sequenceData['n_seq']
                #self.hits = sequenceData['hits']               
                self.hits = []  
                for hit in sequenceData['hits']:
                        newHit = Hit(hit)
                        self.hits.append(newHit)                
                self.profile = Profile(self.pseq)
# add a class for the blast hit results (which contains list of the hits)
class Hit:
        def __init__(self, hitData):
                self.num = hitData['num']
                self.description = hitData['description']
                self.length = hitData['len']
                self.hits = []
                for hit in hitData['hsps']:
                        newHit = Hsps(hit)
                        self.hits.append(newHit)
                #self.hsps = Hsps(hspsData)


# define a class for the high scoring results of blast
## add extra characteristics to this class as/when neccesary
class Hsps:
        def __init__(self, hsps):
                self.dat = hsps


#define profile class to create N-mer profile
class Profile:
    def __init__(self, sequence):
        self.sequence = sequence
        self.kmers = {}
        self.length = len(self.sequence)
        for i in range(0, len(self.sequence)-6+1):
            kmer = self.sequence[i:i+6]
            if kmer not in self.kmers:
                self.kmers[kmer] = 1
            else:
                self.kmers[kmer] += 1
        self.sortedKmers = sorted(self.kmers.items(), reverse=True, key=lambda x: x[1])
        self.entropy = Entropy(self.length, self.kmers)
# define entropy class to calculate evenness
class Entropy:
    def __init__(self, length, kmers):
        totalkmers = length - 6 + 1
        H = 0
        for freq in kmers.values():
            pi = freq/totalkmers
            piLogPi = pi * math.log2(pi)
            H = H + piLogPi
        self.H = - H
        try:
            self.HNorm = self.H / math.log2(totalkmers)
        except (ZeroDivisionError, ValueError):
            self.HNorm = 1

#############################################################################################

# makes a list to add the species objects to
speciesList = []

# try to open each folder in database directory 
for folder in os.listdir(rootDir):
        try:
                with open("%s%s/db.json"%(rootDir,folder)) as file:
                        #create a new object for each species, using json data from the file and add to the list of species
                        data = json.load(file)
                        newSpec = Species(data)
                        speciesList.append(newSpec)
        except NotADirectoryError:
                pass
speciesList.sort(key=lambda x: x.taxon)
numbSpecies = len(speciesList)

def getTotalSequeces():
    numbSeqs = 0
    for spec in speciesList:
        seqs = len(spec.seqs)
        numbSeqs += seqs
    return numbSeqs

#############################################################################################


## Defines a class which will create main widget
class home_Window(QMainWindow):


    # constructor for the class
    def __init__(self):
        super().__init__()  # 

        #creation of GUI will be performed by .initUI method
        self.initUI()

    #define method to build UI window with all its features
    def initUI(self):

        # define some global variables
        self.selectedSpeciesSpecList = None
        self.selectedSpeciesHMMList = None
        self.hmmRowNo = None
        self.selectedHmmSeq = None

        # label to show data of clicked list icon
        self.list_info = QLabel(self)
        self.list_info.move(55,560)
        self.list_info.resize(250,125)
        self.list_info.setTextInteractionFlags(Qt.TextSelectableByMouse)
        
        # search box
        self.species_search_box = QLineEdit(self)
        self.species_search_box.move(50,15)
        self.species_search_box.resize(150,30)
        # search button - filter the species list
        self.filter_btn = QPushButton("Filter", self)
        self.filter_btn.move(210,15)
        self.filter_btn.resize(70,30)
        self.filter_btn.clicked.connect(self.load_spec_list)
        #clear button - clear the list filter
        self.clear_filter_btn = QPushButton("Clear", self)
        self.clear_filter_btn.move(300,15)
        self.clear_filter_btn.resize(70,30)
        self.clear_filter_btn.clicked.connect(self.clear_btn_clicked)
        # Create a list widget to display list of species
        self.spec_list = QListWidget(self)
        self.spec_list.move(50,50)
        self.spec_list.resize(250,500)
        # add each species taxon name to the list
        self.load_spec_list()           
        # call the function to change the text in 'list-Info box' when list item clicked
        self.spec_list.currentItemChanged.connect(self.spec_list_item_changed)
        # make a push button to retrieve species sequences
        self.get_spec_data_btn = QPushButton("Get Data", self)
        self.get_spec_data_btn.move(50,685)
        # action of button to display the data in next list (adds each hmm sequence)
        self.get_spec_data_btn.clicked.connect(self.get_data_clicked)
        # list to display hmm sequences
        self.selectedHmmSeq_list = QListWidget(self)
        self.selectedHmmSeq_list.move(320,110)
        self.selectedHmmSeq_list.resize(500,250)
        self.selectedHmmSeq_title = QLabel(self)
        self.selectedHmmSeq_title.move(320,50)
        self.selectedHmmSeq_title.resize(500,50)
        # Scroll area to print info of each sequences
        self.seq_info_scroll = QScrollArea(self)
        self.seq_info_scroll.move(320, 375)
        self.seq_info_scroll.resize(800,320)
        self.seq_info = QLabel(self)
        self.seq_info.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.seq_info.setWordWrap(True)
        self.seq_info_scroll.setWidget(self.seq_info)
        self.seq_info_scroll.setWidgetResizable(True)
        layout = QVBoxLayout(self.seq_info)
        layout.setAlignment(Qt.AlignTop)
        self.seq_info.setContentsMargins(30,30,30,30)
        self.seq_info.resize(300,320)

        # action to connect on clicked
        self.selectedHmmSeq_list.currentItemChanged.connect(self.hmm_list_item_clicked)

        #button obtain blast hits
        self.get_seq_data_btn = QPushButton("See BLAST Hits", self)
        self.get_seq_data_btn.move(320,700)
        self.get_seq_data_btn.resize(250,25)
        #action when button is pressed
        self.get_seq_data_btn.clicked.connect(self.blast_btn_clicked)

        self.database_info_box = QLabel(self)
        self.database_info_box.move(75,770)
        self.database_info_box.resize(500,40)
        self.database_info_box.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        self.database_info_box.setText("Species:\t%d\t\t\tTotal Sequences:\t%d"%(numbSpecies, getTotalSequeces()))



        #sizing of window 
        self.resize(1250,800)
        self.center() #calls center function as defined below
        self.setWindowTitle("SILKBASE\t\t-\t\tSpider Silk Database")
        # show the window
        self.show()

    
    #function to centre window
    def center(self):        
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())
   
    @pyqtSlot()
    def load_spec_list(self):
        boxText = self.species_search_box.text()
        self.spec_list.clear()
        if boxText == "":
            for spec in speciesList:
                self.spec_list.addItem(spec.taxon)
        else:
            for spec in speciesList:
                if boxText.lower() in spec.taxon.lower() or boxText.lower() in spec.accession.lower():
                    self.spec_list.addItem(spec.taxon)
                else:
                    condition = False
                    for seq in spec.seqs:
                        if boxText.lower() in seq.type:
                            condition = True
                    if condition == True:
                        self.spec_list.addItem(spec.taxon)

    @pyqtSlot()
    def clear_btn_clicked(self):
        self.species_search_box.clear()
        self.load_spec_list()

    # action for clicking species list item
    @pyqtSlot()
    def spec_list_item_changed(self):
        try:
            specName = self.spec_list.currentItem().text()
            for spec in speciesList:
                if specName == spec.taxon:
                    self.change_selected_spec(spec) 
        except:
            pass

    @pyqtSlot()
    def change_selected_spec(self, spec):
        self.selectedSpeciesSpecList = spec
        self.list_info.setText("self.selectedSpeciesSpecList.taxon")
        accession = self.selectedSpeciesSpecList.accession
        taxon = self.selectedSpeciesSpecList.taxon
        noSeqs = len(self.selectedSpeciesSpecList.seqs)
        self.list_info.setText("Taxon:\t%s\n\nAccession:\t%s\n\nSequences:\t%d"%(taxon, accession, noSeqs))

    # action for clicking species list get data button    
    @pyqtSlot()
    def get_data_clicked(self):
        self.selectedSpeciesHMMList = self.selectedSpeciesSpecList
        self.selectedHmmSeq_title.setText(self.selectedSpeciesHMMList.taxon + "\n\nSequences:")
        self.selectedHmmSeq_list.clear()
        for seq in self.selectedSpeciesHMMList.seqs:
                self.selectedHmmSeq_list.addItem(seq.title)
        
    # action for clicking item of hmm list
    @pyqtSlot()
    def hmm_list_item_clicked(self):
        self.hmmRowNo = self.selectedHmmSeq_list.currentRow()
        self.selectedHmmSeq = self.selectedSpeciesHMMList.seqs[self.hmmRowNo]
        self.seq_info.setText("Title:\t%s\n\nType:\t%s\n\nP Seq:\t%s\nLength:\t%s\nEvenness:\t%s\n\nN Seq:\t%s\nLength: %d\n\nNo of BLAST hits:\t%d"
            %(self.selectedHmmSeq.title, self.selectedHmmSeq.type, self.selectedHmmSeq.pseq, len(self.selectedHmmSeq.pseq), self.selectedHmmSeq.profile.entropy.HNorm, self.selectedHmmSeq.nseq, len(self.selectedHmmSeq.nseq),len(self.selectedHmmSeq.hits)))
        

    @pyqtSlot()
    def blast_btn_clicked(self):
        self.blast_result_window = BLASTPopup()
        self.selectedHmmSeq = self.selectedSpeciesHMMList.seqs[self.hmmRowNo]
        self.blast_result_window.setWindowTitle("BLAST Hits")
        printStr = ""
        for Hit in self.selectedHmmSeq.hits:
            for desc in Hit.description:
                printStr = printStr + str(desc).replace("'","").replace(",","\n\n").replace(":",":\t\t").replace("{","").replace("}","") + "\n\n"
            printStr = printStr + (
                 "Length:\t\t" + str(Hit.length) +"\n\n"
                + "Fragment Number:\t\t" + str(Hit.num) + "\n\n")
            printStr = printStr + "Number of hits:\t" + str(len(Hit.hits)) + "\n\n"
            for hsps in Hit.hits:
                printStr = printStr + "\t" + str(hsps.dat).replace("'","").replace(",","\n\t").replace(":",":\t\t").replace("{","").replace("}","") + "\n\n\n\n"

        self.blast_result_window.textbox.setText(printStr)
        self.blast_result_window.show()

class BLASTPopup(QScrollArea):
    def __init__(self):
        super(BLASTPopup, self).__init__()
        self.textbox = QLabel(self)
        self.textbox.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.setWidgetResizable(True)
        layout = QVBoxLayout(self.textbox)
        layout.setAlignment(Qt.AlignTop)
        self.setWidget(self.textbox)
        self.resize(900,700)
        self.show()


# instantiates window 
if __name__  == '__main__':
    
    #create a new QApplication 
    app = QApplication(sys.argv)
    
    wind = home_Window()
    sys.exit(app.exec_())

