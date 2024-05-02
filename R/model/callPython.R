source(file='common.R')
print(getwd())

library(reticulate)


env <- globalenv()

#install_miniconda(path = "C:\\Users\\magnu\\miniconda3", update = T)
use_condaenv(condaenv = "dml", conda = "C:/Users/magnu/miniconda3/condabin/conda.bat")
#use_condaenv(condaenv = "r-reticulate", conda = "C:\\Users\\magnu\\python\\envs\\r-reticulate/python.exe")


#pd <- import("pandas")


cwd <- getwd()

library(os)

setwd("C:\\Users\\magnu\\OneDrive - Västra Götalandsregionen\\DIRE\\model\\dire")
source_python("C:\\Users\\magnu\\OneDrive - Västra Götalandsregionen\\DIRE\\model\\bin\\imports.py")

BERTModel()
setwd(cwd)

vocab <- torch.load(vocabfile_name)

nn_net = BERTModel(vocab_size=len(vocab),
                   num_hiddens=config['parameters']['network']['num_hiddens'],
                   norm_shape=config['parameters']['network']['norm_shape'],
                   ffn_num_input=config['parameters']['network']['ffn_num_input'],
                   ffn_num_hiddens=config['parameters']['network']['ffn_num_hiddens'],
                   attention_heads=config['parameters']['network']['attention_heads'],
                   attention_layers=config['parameters']['network']['attention_layers'],
                   dropout=config['parameters']['network']['dropout'],
                   key_size=config['parameters']['network']['key_size'],
                   query_size=config['parameters']['network']['query_size'],
                   value_size=config['parameters']['network']['value_size'],
                   number_ab=len(config['antibiotics']))



#data = pd.DataFrame({'Antibiotic':
#    ['AMP_S,AMX_S,AMC_S,PIP_S,TZP_S,CAZ_S,CRO_S,CTX_S,FEP_S,CIP_S,OFX_S,LVX_S,MFX_S,GEN_S,TOB_S',
#     'AMP_R,AMX_S,AMC_S,PIP_R,TZP_S,CAZ_S,CRO_S,CTX_S,FEP_R,CIP_S,OFX_S,LVX_S,MFX_S,GEN_S,TOB_S'],
#  'n': [15, 15], 'x': ['SE 72 F 2018-10', 'SE 63 M 2016-03']})






