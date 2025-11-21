print("Find data and vocab")      
print(config)

vocabfile_name = os.path.join(DIRE_PATH, config['vocab']['vocab'])
print(vocabfile_name)  

print("Start process")  
vocab = torch.load(vocabfile_name)
vocab_size = len(vocab.idx_to_token)
config["vocab"]["vocab_size"] = vocab_size    
PAD_id = vocab.token_to_idx['<pad>']
Fake_position = config["vocab"]["pad_position"] # pad position
pos_response = config["model"]["pos_response"]  # 1 for resistant
ab2 = list(pos_response.keys())                 # the list of possible antibiotics
ab2 = [l.split("_")[0] for l in ab2]            # the list of possible antibiotics
n_metadata = config["train"]["n_metadata"] # (so cls does count)
n_antibiotics = config["model"]["pos_response"]  
n_antibiotics = len(list(n_antibiotics.keys()))  

model_folder = "model_with_patient_data"
data_model = "both"

print("Started process")      

def INIT():
    print("Start INIT")  
    #combination=('AMC', 'TZP', 'PIP', 'CRO', 'CTX')
    #print(pCombinations)
    rank = 0
    data_models = data_model
    model_antibiotics =  create_transformer_model(config, PAD_id, Fake_position, rank, data = "antibiotics")
    if data_models == "both":
        model_patient_info =  create_transformer_model(config, PAD_id, Fake_position, rank, data = "patient")
        classification_model = create_classification_model(config, rank, "both")
    else:
        classification_model = create_classification_model(config, rank, "antibiotics")
        model_patient_info = None
    #
    device_id = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    state = load_previous_model(model_folder, model_antibiotics, model_patient_info,
                            classification_model, rank, load_previous=True,
                            data=data_models, device=device_id)
    model_antibiotics.to(device_id)
    if data_models == "both":
        model_patient_info.to(device_id)
    classification_model.to(device_id)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")


    print("Finished INIT") 
    return data_models,model_antibiotics, model_patient_info, classification_model, device, rank

def toInputFile(pData,pCombinations):
    data = pData
    combinations = pCombinations
    combos = combinations.tolist()
   # Split antibiotics for each row
    rows = []
    for _, row in data.iterrows():
     abx = [a.strip() for a in row["Antibiotic"].split(",")]
     combo_group = [a for a in abx if any(a.startswith(c) for c in combos)]
     other_group = [a for a in abx if a not in combo_group]
     rows.append(["ESCCOL", row["x"], " ".join(combo_group), " ".join(other_group)])

   # Create final DataFrame (no headers)
    input_df = pd.DataFrame(rows) 
   

    #print(input_df)
    
#    data_combo.to_csv("data_combo.csv",FALSE)
#    data_other.to_csv("data_other.csv",FALSE)

    
    #inputData = pd.DataFrame(rep()
    return input_df
    

 

def RUN(pData,pSignificant_levels,pCombinations,pData_models,pModel_antibiotics, pModel_patient_info, pClassification_model, pDevice, pRank):
    print("Start RUN")  
#    assert (0. <= pSignificant_levels <= 1.)
    data = pData
#    print(pSignificant_levels)
    significant_levels = pSignificant_levels
#        data = pd.DataFrame({'Antibiotic':
#                                 ['AMP_S,AMX_S,AMC_S,PIP_S,TZP_S,CAZ_S,CRO_S,CTX_S,FEP_S,CIP_S,OFX_S,LVX_S,MFX_S,GEN_S,TOB_S',
#                                  'AMP_R,AMX_S,AMC_S,PIP_R,TZP_S,CAZ_S,CRO_S,CTX_S,FEP_R,CIP_S,OFX_S,LVX_S,MFX_S,GEN_S,TOB_S'],
#                             'n': [15, 15], 'x': ['SE 72 F 2018-10', 'SE 63 M 2016-03']})

#    data_val = AntibioticDataset(df=data, encode_response=encode_response, vocab=vocab, combination=combinations,
#                                 x_max_len=config['parameters']['sampler']['x_max_len'],
#                                 y_max_len=config['parameters']['sampler']['y_max_len'])
#    data_val = DataLoader(data_val, batch_size=min(config['parameters']['batch_size'], data.shape[0]),
#                          shuffle=False)
#
#    given, preds, answer, softmax, abs_, pred_index = model.compile_results(data_val)
#    given = given.astype(str).replace('nan', 'predict').replace('0.0', 'S').replace('1.0', 'R')
#    preds = preds.astype(str).replace('nan', 'given').replace('0.0', 'S').replace('1.0', 'R')
#    answer = answer.astype(str).replace('nan', 'given').replace('0.0', 'S').replace('1.0', 'R')
#    conformal_preds = list()
#    for level in significant_levels:
#      cp = ConformalPredictionAntibiotic(antibiotics=config['antibiotics'],alphas=pickle.load(open("saved_models/{}.p".format(config['load']['conformal_name']),"rb")),significant_level=level)
#      conformal_pred = cp.compile_results(softmax, abs_, pred_index)
#      conformal_preds.append(conformal_pred)

#    print(pData)
#    print(pCombinations)
    inputFile = toInputFile(pData,pCombinations)
    csv_buffer = StringIO()
    inputFile.to_csv(csv_buffer, index=False, header=False)
    csv_buffer.seek(0)
    sentence = read_data(csv_buffer, n_metadata, n_antibiotics)
    
#    print(sentence)
    
    df_test_f = forward_val(sentence, config, pModel_antibiotics, pModel_patient_info, pClassification_model, pDevice, pRank, 
    vocab, Fake_position, pData_models)
    
    #print(df_test_f)
    
    print("Finished RUN")  
    return df_test_f
