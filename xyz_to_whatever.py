
from tkinter import *
from tkinter import ttk as ttk



def _Preview(*args, **kwargs):
    global Preview, Method, CalcMain, CC, CCSOPPA, CalcAdd, CalcRes
    Preview = '''**DALTON INPUT'''
    try:
        CalcMain
    except NameError:
        if CalcRes == '':
            Preview += '''\n.RUN WAVEFUNCTIONS'''
        else: 
            Preview += '''\n.RUN RESPONSE'''
    else:
        if CalcMain == '.RUN WAVEFUNCTIONS' and CalcRes != '':
            Preview += '''\n.RUN RESPONSE'''
        else:
            Preview += f'''\n{CalcMain}'''

    try:
        Method
    except NameError:
        Preview += '''\n**WAVEFUNCTIONS\n.HF'''
    else:
        Preview += '''\n**WAVEFUNCTIONS'''
        if Method == 'HF':
            Preview += '''\n.HF'''
        elif Method == 'MP2':
            Preview += '''\n.HF\n.MP2'''
        elif Method == 'CC':
            Preview += '''\n.HF\n.CC'''
        elif Method == 'DFT':
            try:
                Functional
            except NameError:
                Preview += f'''\n.DFT\nLDA'''
            else:
                Preview += f'''\n.DFT\n{Functional}'''
        else:
            Preview += '''\n.HF'''
        
        if Method == 'CC':
            Preview += '''\n*CC INPUT'''
            for i in CC:
                Preview += f'''\n{i}'''
            try:
                CCSOPPA
            except NameError:
                pass
            else:
                if CCSOPPA != '':
                    Preview += f'''\n.{CCSOPPA}'''
    
    if '**PROPERTIES' in CalcAdd:
        Preview += '''\n**PROPERTIES'''

    if '**RESPONSE' in CalcAdd:
        Preview += '''\n**RESPONSE'''
    
    Preview += '''\n*END OF INPUT'''
    
    _PreviewText()

def _PreviewText(*args, **kwargs):
    global PreviewLabel
    PreviewLabel.config(font=('Times New Roman', 12), text = Preview, justify='left')
    PreviewLabel.grid(row=0, column=0)

def _CalcType(*args, **kwargs):
    global CalcMain, CalcTypeMain
    CalcMain = CalcTypeMain.get()
    _Preview()

def _CalcTypeAdd(*args, **kwargs):
    global CalcTypeProp, CalcTypeRes, CalcAdd, tabControl, CalcProp, CalcRes
    CalcProp = CalcTypeProp.get()
    CalcRes = CalcTypeRes.get()
    CalcAdd = [i.get() for i in [CalcTypeProp, CalcTypeRes] if i.get() != '']
    _TabRemove()

    if CalcProp != '':
        tabControl.tab(2, state='normal')
    if CalcRes != '':
        tabControl.tab(3, state='normal')

    _Preview()

def _TabRemove(*args, **kwargs):
    global tabControl
    if CalcProp == '':
        tabControl.tab(2, state='hidden')
    if CalcRes == '':
        tabControl.tab(3, state='hidden')

def _MethodType(*args, **kwargs):
    global Method
    Method = MethodType.get()
    _MethodRemove()

    # Creation of Functional type dropdown box
    if Method == 'DFT':
        FunctionalLabel.grid(row=1, column=3, padx=10, pady=10)
        FunctionalTypes.grid(row=1 ,column=4)
        return
    if Method == 'CC':
        CCLabel.grid(row=2, column=1, padx=3, pady=10, sticky=W)
        CCSOPPALabel.grid(row=3, column=1, padx=3, pady=10, sticky=W)
        CCTypesCCS.grid(row=2, column=2, padx=3, pady=10, sticky=W)
        CCTypesMP2.grid(row=2, column=3, padx=3, pady=10, sticky=W)
        CCTypesCC2.grid(row=2, column=4, padx=3, pady=10, sticky=W)
        CCTypesCISpD.grid(row=2, column=5, padx=3, pady=10, sticky=W)
        CCTypesCCSD.grid(row=2, column=6, padx=3, pady=10, sticky=W)
        CCTypesCCSDRp3.grid(row=2, column=7, padx=3, pady=10, sticky=W)
        CCTypesCCSDpT.grid(row=2, column=8, padx=3, pady=10, sticky=W)
        CCTypesCC3.grid(row=2, column=9, padx=3, pady=10, sticky=W)
        CCTypesSOPPA.grid(row=3, column=2, padx=3, pady=10, columnspan=3)
        return

    _Preview()

def _MethodRemove(*args, **kwargs):
    if Method != 'DFT':
        FunctionalTypes.grid_remove()
        FunctionalLabel.grid_remove()

    if Method != 'CC':
        CCLabel.grid_remove()
        CCSOPPALabel.grid_remove()
        CCTypesCCS.grid_remove()
        CCTypesMP2.grid_remove()
        CCTypesCC2.grid_remove()
        CCTypesCISpD.grid_remove()
        CCTypesCCSD.grid_remove()
        CCTypesCCSDRp3.grid_remove()
        CCTypesCCSDpT.grid_remove()
        CCTypesCC3.grid_remove()
        CCTypesSOPPA.grid_remove()

def _FunctionalType(*args, **kwargs):
    global Functional
    if Method == 'DFT':
        Functional = FuncType.get()
        _Preview()
        return

# def _BasissetStyles(*args, **kwargs):
#     global BasissetStyle
#     BasissetStyle = BasissetStyles.get()

def _CCType(*args, **kwargs):
    global CC, CCSOPPA
    if Method == 'CC':
        CC = [i.get() for i in [CCTypeCCS, CCTypeMP2, CCTypeCC2, CCTypeCISpD, CCTypeCCSD, CCTypeCCSDRp3, CCTypeCCSDpT, CCTypeCC3] if i.get() != '']
        CCSOPPA = CCTypeSOPPA.get()
        _Preview()
        return
    pass

if __name__ == '__main__':
    root = Tk()
    root.title('XYZ to DALTON20')
    tabControl = ttk.Notebook(root)

    root.geometry('800x500')

    # Tabs can be added here
    MainTab = Frame(tabControl)
    MethodTab = Frame(tabControl)
    PreviewTab = Frame(tabControl)

    PropertiesTab = Frame(tabControl)
    ResponseTab = Frame(tabControl)

    tabControl.add(MainTab, text='Main')
    tabControl.add(MethodTab, text='Method')
    tabControl.add(PropertiesTab, text='Properties', state='hidden')
    tabControl.add(ResponseTab, text='Response', state='hidden')
    tabControl.add(PreviewTab, text='Preview')
    tabControl.pack(expand=1, fill="both")

    # Creation of Preview
    Preview = '''**DALTON INPUT'''
    PreviewLabel = Label(PreviewTab)
    PreviewLabel.config(font=('Times New Roman', 12),text = Preview)
    PreviewLabel.pack()
    # _Preview()

    # Creation of Functional type dropdown box
    FunctionalLabel = Label(MethodTab, text='Functional', font = ("Times New Roman", 12))

    FuncType = StringVar()

    FunctionalTypes = ttk.Combobox(MethodTab, textvariable=FuncType, values=('LDA', 'BLYP', 'B3LYP','CAMB3LYP', 'B2PLYP', 'PBE', 'PBE0', 'PBE0DH'), state='readonly')
    FunctionalTypes.set('LDA')
    FunctionalTypes.bind("<<ComboboxSelected>>", _FunctionalType)
    
    # Creation of CC type radiobutton and dropdown boxex
    CCLabel = Label(MethodTab, text='CC module inputs:', font = ("Times New Roman", 12))
    CCSOPPALabel = Label(MethodTab, text='CC SOPPA inputs', font = ("Times New Roman", 12))

    CCTypeCCS = StringVar()
    CCTypeMP2 = StringVar()
    CCTypeCC2 = StringVar()
    CCTypeCISpD = StringVar()
    CCTypeCCSD = StringVar()
    CCTypeCCSDRp3 = StringVar()
    CCTypeCCSDpT = StringVar()
    CCTypeCC3 = StringVar()

    CCTypeSOPPA = StringVar()
    
    CCTypesCCS = ttk.Checkbutton(MethodTab, text='CCS', onvalue='.CCS', offvalue='', var=CCTypeCCS, command=_CCType)
    CCTypesMP2 = ttk.Checkbutton(MethodTab, text='MP2', onvalue='.MP2', offvalue='', var=CCTypeMP2, command=_CCType)
    CCTypesCC2 = ttk.Checkbutton(MethodTab, text='CC2', onvalue='.CC2', offvalue='', var=CCTypeCC2, command=_CCType)
    CCTypesCISpD = ttk.Checkbutton(MethodTab, text='CCS(D)', onvalue='.CCS(D)', offvalue='', var=CCTypeCISpD, command=_CCType)
    CCTypesCCSD = ttk.Checkbutton(MethodTab, text='CCSD', onvalue='.CCSD', offvalue='', var=CCTypeCCSD, command=_CCType)
    CCTypesCCSDRp3 = ttk.Checkbutton(MethodTab, text='CCSDR(3)', onvalue='.CCSDR(3)', offvalue='', var=CCTypeCCSDRp3, command=_CCType)
    CCTypesCCSDpT = ttk.Checkbutton(MethodTab, text='CCSD(T)', onvalue='.CCS(T)', offvalue='', var=CCTypeCCSDpT, command=_CCType)
    CCTypesCC3 = ttk.Checkbutton(MethodTab, text='CC3', onvalue='.CC3', offvalue='', var=CCTypeCC3, command=_CCType)

    CCTypesSOPPA = ttk.Combobox(MethodTab, textvariable=CCTypeSOPPA, values=['','SOPPA', 'SOPPA2', 'SOPPA(CCSD)', 'AO-SOPPA'], state='readonly')
    CCTypesSOPPA.bind("<<ComboboxSelected>>", _CCType)

    CC = []

    # Creation of Job type radiobuttons
    CalcLabel = Label(MainTab, text='Calculation Type', font = ("Times New Roman", 12))
    CalcLabel.grid(row=1, column=0, padx=10, pady=10, sticky=W)

    CalcTypeMain = StringVar()

    CalcTypesSingle = ttk.Radiobutton(MainTab, variable=CalcTypeMain, value='.RUN WAVEFUNCTION', text='Single Point', command=_CalcType)
    CalcTypesOptimize = ttk.Radiobutton(MainTab, variable=CalcTypeMain, value='.OPTIMIZE', text='Optimize', command=_CalcType)
    CalcTypeMain.set('.RUN WAVEFUNCTION')

    CalcTypesSingle.grid(row=1, column=1, padx=10, pady=10)
    CalcTypesOptimize.grid(row=1, column=2, padx=10, pady=10)

    # Creation of additional Job type checkbuttons
    CalcLabelAdd = Label(MainTab, text='Additional calculations', font = ("Times New Roman", 12))
    CalcLabelAdd.grid(row=2, column=0, padx=10, pady=10)

    CalcTypeProp = StringVar()
    CalcTypeRes = StringVar()

    CalcTypesProp = ttk.Checkbutton(MainTab, text='Properties', onvalue='**PROPERTIES', offvalue='', var=CalcTypeProp, command=_CalcTypeAdd)
    CalcTypesRes = ttk.Checkbutton(MainTab, text='Response', onvalue='**RESPONSE', offvalue='', var=CalcTypeRes, command=_CalcTypeAdd)

    CalcTypesProp.grid(row=2, column=1, padx=10, pady=10, sticky=W)
    CalcTypesRes.grid(row=2, column=2, padx=10, pady=10, sticky=W)

    CalcAdd = []
    CalcRes = ''
    CalcProp = ''
    
    # Creation of Method type dropdown box
    MethodLabel = Label(MethodTab, text='Method', font = ("Times New Roman", 12))
    MethodLabel.grid(column=0, row=1, padx=10, pady=10)

    MethodType = StringVar()

    MethodTypes = ttk.Combobox(MethodTab, textvariable=MethodType, values=('HF', 'MP2', 'CC', 'DFT'), state='readonly')
    MethodTypes.set('HF')
    MethodTypes.grid(row=1, column=1, columnspan=2)
    MethodTypes.bind("<<ComboboxSelected>>", _MethodType)

    # # Creation of Basis set type dropdown boxes
    # global BasissetType, BasissetLabel, BasissetStyles, BasissetStyles_dict
    # BasissetStyles_dict = {
    # 'STO':
    #     ('2G', '3G', '6G'),
    # 'Pople': {
    #     '3-21G':
    #         (('',('','*')),
    #         ('++',('','*'))),
    #     '4-31G':
    #         (('','')),
    #     '6-31G':
    #         (('',('','*','**','3df,3pd')),
    #         ('+',('','*')),
    #         ('++',('','*','**'))),
    #     '6-311G':
    #         (('',('','*','**', '(2df,2pd)')),
    #         ('+',('*')),
    #         ('++',('**','(2d,2p)','(3df,3pd)')))
    #     },
    # 'ano':
    #     ('1','2','3','4'),
    # 'Dunniger': {
    #     'cc-pVDZ': ('','aug'),
    #     'cc-pVTZ': ('','aug'),
    #     'cc-pVQZ': ('','aug'),
    #     'cc-pV5Z': ('','aug'),
    #     'cc-pV6Z': ('','aug')
    #     },
    # 'pc': {
    #     '1': ('','aug'),
    #     '2': ('','aug'),
    #     '3': ('','aug'),
    #     '4': ('','aug')
    #     }
    # }
    # BasissetLabel = Label(MethodTab, text='Basis set type', font = ("Times New Roman", 12))
    # BasissetLabel.grid(row=3, column=0, padx=10, pady=10)

    # BasissetType = StringVar()

    # BasissetStyles = ttk.Combobox(MethodTab, textvariable=BasissetType, values=('STO', 'Pople', 'ano', 'Dunninger', 'Polarization consistent', ), state='readonly')
    # BasissetStyles.set('STO')
    # BasissetStyles.grid(row=3, column=1)
    # BasissetStyles.bind("<<ComboboxSelected>>", _BasissetStyles)


    _Preview()

    print(Preview)