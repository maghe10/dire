data_folder = os.path.join(ROOT_DIR, 'data')
vocabfile_name = os.path.join(data_folder, config['files']['vocab'])

## FIXME make class
encode_response = encoding_resp(config['antibiotics'], config['encode_respone'])
vocab = torch.load(vocabfile_name)
#model = 0
#combination = 0

def INIT(pCombinations):
    #combination=('AMC', 'TZP', 'PIP', 'CRO', 'CTX')
    combinations = pCombinations
    devices = d2l.try_all_gpus()
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

    nn_net.to(devices[0])
    nn_net.to_within(devices[0])

    model = AntibioticModelEval(net=nn_net,
                                devices=devices,
                                x_interval=[len(combinations), len(combinations)],
                                antibiotics=config['antibiotics'])
    model.load_model(config['load']['model_name'])
    return model, combinations

def RUN(pData,pSignificant_levels,pModel,pCombinations):
#    assert (0. <= pSignificant_levels <= 1.)
    data = pData
#    print(pSignificant_levels)
    significant_levels = pSignificant_levels
    model = pModel
    combinations = pCombinations
#        data = pd.DataFrame({'Antibiotic':
#                                 ['AMP_S,AMX_S,AMC_S,PIP_S,TZP_S,CAZ_S,CRO_S,CTX_S,FEP_S,CIP_S,OFX_S,LVX_S,MFX_S,GEN_S,TOB_S',
#                                  'AMP_R,AMX_S,AMC_S,PIP_R,TZP_S,CAZ_S,CRO_S,CTX_S,FEP_R,CIP_S,OFX_S,LVX_S,MFX_S,GEN_S,TOB_S'],
#                             'n': [15, 15], 'x': ['SE 72 F 2018-10', 'SE 63 M 2016-03']})

    data_val = AntibioticDataset(df=data, encode_response=encode_response, vocab=vocab, combination=combinations,
                                 x_max_len=config['parameters']['sampler']['x_max_len'],
                                 y_max_len=config['parameters']['sampler']['y_max_len'])
    data_val = DataLoader(data_val, batch_size=min(config['parameters']['batch_size'], data.shape[0]),
                          shuffle=False)

    given, preds, answer, softmax, abs_, pred_index = model.compile_results(data_val)
    given = given.astype(str).replace('nan', 'predict').replace('0.0', 'S').replace('1.0', 'R')
    preds = preds.astype(str).replace('nan', 'given').replace('0.0', 'S').replace('1.0', 'R')
    answer = answer.astype(str).replace('nan', 'given').replace('0.0', 'S').replace('1.0', 'R')
    conformal_preds = list()
    for level in significant_levels:
      cp = ConformalPredictionAntibiotic(antibiotics=config['antibiotics'],alphas=pickle.load(open("saved_models/{}.p".format(config['load']['conformal_name']),"rb")),significant_level=level)
      conformal_pred = cp.compile_results(softmax, abs_, pred_index)
      conformal_preds.append(conformal_pred)
      
    return given, preds, answer, conformal_preds
