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
    return input_df
    

def RUN_INPUT_FRAME(pInputFrame,pData_models,pModel_antibiotics, pModel_patient_info, pClassification_model, pDevice, pRank):
    #if isinstance(pData, pd.DataFrame) and len(pData) > 0:
    #            first_row_str = pData.iloc[[0]].to_string(index=False)
    #else:
    #            first_row_str = "<no DataFrame or empty>"

    #print(
    #    f"""
    #    === RUN DEBUG =====================================
    #    Combinations : {list(pCombinations)}
    #    Significance : {pSignificant_levels}
    #    Device       : {pDevice}
    #    Model-type   : {pData_models}
    #    Rank         : {pRank}
    #    Data rows    : {len(pData)}
    #    First row:
    #    {first_row_str}
    #    ===================================================================
    #    """
    #)
    print("Start RUN")  

    inputFile = pInputFrame
    #print(inputFile)
    #Stop()
    csv_buffer = StringIO()
    inputFile.to_csv(csv_buffer, index=False, header=False)
    csv_buffer.seek(0)
    sentence = read_data(csv_buffer, n_metadata, n_antibiotics)
    #print(sentence)
    df_test_f = forward_val(sentence, config, pModel_antibiotics, pModel_patient_info, pClassification_model, pDevice, pRank, 
    vocab, Fake_position, pData_models)
    #print(df_test_f)
    print("Finished RUN")  
    return df_test_f
 

def RUN(pData,pSignificant_levels,pCombinations,pData_models,pModel_antibiotics, pModel_patient_info, pClassification_model, pDevice, pRank):
    #if isinstance(pData, pd.DataFrame) and len(pData) > 0:
    #            first_row_str = pData.iloc[[0]].to_string(index=False)
    #else:
    #            first_row_str = "<no DataFrame or empty>"

    #print(
    #    f"""
    #    === RUN DEBUG =====================================
    #    Combinations : {list(pCombinations)}
    #    Significance : {pSignificant_levels}
    #    Device       : {pDevice}
    #    Model-type   : {pData_models}
    #    Rank         : {pRank}
    #    Data rows    : {len(pData)}
    #    First row:
    #    {first_row_str}
    #    ===================================================================
    #    """
    #)
    print("Start RUN")  
    data = pData
    significant_levels = pSignificant_levels

    inputFile = toInputFile(pData,pCombinations)
    #print(inputFile)
    #Stop()
    csv_buffer = StringIO()
    inputFile.to_csv(csv_buffer, index=False, header=False)
    csv_buffer.seek(0)
    sentence = read_data(csv_buffer, n_metadata, n_antibiotics)
    #print(sentence)
    df_test_f = forward_val(sentence, config, pModel_antibiotics, pModel_patient_info, pClassification_model, pDevice, pRank, 
    vocab, Fake_position, pData_models)
    #print(df_test_f)
    print("Finished RUN")  
    return df_test_f
