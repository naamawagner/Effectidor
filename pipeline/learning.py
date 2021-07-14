import pandas as pd
import numpy as np
import csv
import matplotlib.pyplot as plt
from sklearn import model_selection
from sklearn.metrics import average_precision_score
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.feature_selection import SelectFromModel, RFECV
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import os
import seaborn as sns
#%%
out_dir='out_learning'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
features_file=r'features.csv'
dataset=pd.read_csv(features_file)
dataset.is_effector.replace(['no','effector'],[0,1],inplace=True)
array = dataset.values
features = array[:,1:-1]
un_labeled=dataset[dataset.is_effector=='?']
labeled=dataset[dataset.is_effector!='?']

labeled.is_effector.astype('int')
yes= labeled[labeled.is_effector==1]
no=labeled[labeled.is_effector==0]

amount_of_positives = len(yes)
K=10 # number of folds
if 0.8*amount_of_positives < 10:
    K = min(int(0.8*amount_of_positives),5)
seed = 7
#%%
all_genome=dataset.values
X_all=all_genome[:,1:-1]
Y_all=all_genome[:,-1]
y=pd.DataFrame(Y_all)
names=all_genome[:,0]
genes=pd.DataFrame(names)

#normalization
norm_scaler = MinMaxScaler(feature_range=(0, 1)) 
normalizedALL = norm_scaler.fit_transform(X_all)
norm=pd.DataFrame(normalizedALL)
norm_data=pd.concat([genes,norm,y],axis=1)
norm_data.columns=dataset.columns
norm_array=norm_data.values
norm_features=norm_array[:,1:-1]
labeled_norm=norm_data[norm_data.is_effector!='?']
array_norm=labeled_norm.values
X_norm=array_norm[:,1:-1]
Y_norm=array_norm[:,-1]
Y_norm=Y_norm.astype('int')
X_train_norm, X_test_norm, y_train_norm, y_test_norm = model_selection.train_test_split(X_norm,Y_norm,test_size=0.2, random_state=seed,stratify=Y_norm)


#standardization
standard_scaler = StandardScaler().fit(X_all)
standardizedALL = standard_scaler.transform(X_all)
stand=pd.DataFrame(standardizedALL)
stand_data=pd.concat([genes,stand,y],axis=1)
stand_data.columns=dataset.columns
stand_array=stand_data.values
stand_features=stand_array[:,1:-1]
labeled_stand=stand_data[stand_data.is_effector!='?']
array_stand=labeled_stand.values
X_stand=array_stand[:,1:-1]
Y_stand=array_stand[:,-1]
Y_stand=Y_stand.astype('int')
X_train_stand, X_test_stand, y_train_stand, y_test_stand = model_selection.train_test_split(X_stand,Y_stand,test_size=0.2, random_state=seed,stratify=Y_stand)


#untouched
array1=labeled.values
X=array1[:,1:-1]
Y=array1[:,-1]
Y=Y.astype('int')
X_train, X_test, y_train, y_test = model_selection.train_test_split(X,Y,test_size=0.2, random_state=seed,stratify=Y)

#%%
scoring = 'average_precision'

models = []
models.append(('LDA', LinearDiscriminantAnalysis()))
models.append(('NB',GaussianNB()))
models.append(('KNN', KNeighborsClassifier(weights='distance')))
models.append(('SVM', SVC(probability=True)))

x_y_options={'untouched':(X_train,X_test,y_train,y_test,features,X,Y),'normalized':(X_train_norm, X_test_norm,y_train_norm, y_test_norm,norm_features,X_norm,Y_norm),'standardized':(X_train_stand, X_test_stand,y_train_stand, y_test_stand,stand_features,X_stand,Y_stand)}
new_models=[]
for name,model in models:
    results = []
    names = []
    all_results={}
    for option in x_y_options:
        X_train,Y_train=x_y_options[option][0],x_y_options[option][2]
        scoring = scoring
        kfold = model_selection.StratifiedKFold(n_splits=K,shuffle=True, random_state=seed)
        cv_results = model_selection.cross_val_score\
        (model, X_train, y=Y_train, cv=kfold, scoring=scoring)
        results.append(cv_results)
        names.append(option)
        all_results[option]=cv_results.mean()

    best_mean=max(all_results.values())
    for key in all_results:
        if all_results[key]==best_mean:
            best=key
            break
    X_train,X_test,Y_train,Y_test,full_features,full_X,full_Y=x_y_options[best]
    new_models.append((name,model,best,X_train,X_test,Y_train,Y_test,full_features,full_X,full_Y))
new_models.append(('LR', LogisticRegression(),'normalized',X_train_norm,X_test_norm,y_train_norm,y_test_norm,norm_features,X_norm,Y_norm)) # LR requires normalized features for convergance
#new_models.append(('RDF', RandomForestClassifier(n_estimators=1000,max_depth=10),'untouched',X_train,X_test,y_train,y_test,features,X,Y)) # RDF can handle different "untouched" data
#%%
#feature selection
updated_models=[]    
for name,model,best,X_TRAIN,X_TEST,Y_TRAIN,Y_TEST,full_features,full_X,full_Y in new_models:
    train_features=[]
    test_features=[]
    full_features_l = [] # featues of all labeled data
    results=[]
    all_data_features=[] # features of all data - both labeled and unlabeled
    try:
        #RecursiveFeatureElimination
        kfold = model_selection.StratifiedKFold(n_splits=K,shuffle=True, random_state=seed)
        selector=RFECV(model, step=0.1, cv=kfold,scoring=scoring)
        selector = selector.fit(X_TRAIN, Y_TRAIN)
        X_RFECV_train = selector.transform(X_TRAIN)
        X_RFECV_test = selector.transform(X_TEST)
        X_RFECV = selector.transform(full_X)
        features_RFECV = selector.transform(full_features)
        train_features.append(X_RFECV_train)
        test_features.append(X_RFECV_test)
        full_features_l.append(X_RFECV)
        all_data_features.append(features_RFECV)
        #SelectFromModel
        fit_model=model.fit(X_TRAIN,Y_TRAIN)
        fit_model.feature_importances_ 
        selectmodel = SelectFromModel(fit_model,prefit=True)
        X_from_model_train=selectmodel.transform(X_TRAIN)
        X_from_model_test = selectmodel.transform(X_TEST)
        X_from_model = selectmodel.transform(full_X)
        features_model = selectmodel.transform(full_features)
        train_features.append(X_from_model_train)
        test_features.append(X_from_model_test)
        full_features_l.append(X_from_model)
        all_data_features.append(features_model)
        #SelectFromExtraTree
        clf = ExtraTreesClassifier(n_estimators=50)
        fit_model=clf.fit(X_TRAIN,Y_TRAIN)
        fit_model.feature_importances_
        selectmodel = SelectFromModel(fit_model,prefit=True)
        X_RDF_train=selectmodel.transform(X_TRAIN)
        X_RDF_test=selectmodel.transform(X_TEST)
        X_RDF = selectmodel.transform(full_X)
        features_RDF = selectmodel.transform(features)
        train_features.append(X_RDF_train)
        test_features.append(X_RDF_test)
        full_features_l.append(X_RDF)
        all_data_features.append(features_RDF)
    except: #no feature importance or coef_
        clf = ExtraTreesClassifier(n_estimators=50)
        fit_model=clf.fit(X_TRAIN,Y_TRAIN)
        fit_model.feature_importances_
        selectmodel = SelectFromModel(fit_model,prefit=True)
        X_RDF_train=selectmodel.transform(X_TRAIN)
        X_RDF_test=selectmodel.transform(X_TEST)
        X_RDF = selectmodel.transform(full_X)
        features_RDF = selectmodel.transform(features)
        train_features.append(X_RDF_train)
        test_features.append(X_RDF_test)
        full_features_l.append(X_RDF)
        all_data_features.append(features_RDF)
    train_features.append(X_TRAIN)
    test_features.append(X_TEST)
    full_features_l.append(full_X)
    all_data_features.append(features)
    for feature in train_features:
        kfold = model_selection.StratifiedKFold(n_splits=K,shuffle=True, random_state=seed)
        cv_results = model_selection.cross_val_score\
        (model, X_TRAIN, y=Y_TRAIN, cv=kfold, scoring=scoring)
        results.append(cv_results.mean())
    best_f=max(results)
    i=results.index(best_f)
    best_features_train=train_features[i]
    best_features_test=test_features[i]
    best_features = full_features_l[i]
    all_data_best=all_data_features[i]
    updated_models.append((name,model,best_features_train,best_features_test,Y_TRAIN,Y_TEST,best_features,all_data_best))
updated_models.append(('RDF', RandomForestClassifier(n_estimators=500,max_depth=10,random_state=seed),X_train,X_test,Y_train,Y_test,X,features))
#%%
#comparing the algorithms after their preprocessing
results_train={}
results_test={}
for name,model,X_TRAIN,X_TEST,Y_TRAIN,Y_TEST,all_X,all_data in updated_models:
    predictions_train_cv=model_selection.cross_val_predict(model,X_TRAIN,Y_TRAIN,cv=kfold,method="predict_proba")
    AUPRC_train_cv=average_precision_score(Y_TRAIN,predictions_train_cv[:,1])
    model.fit(X_TRAIN,Y_TRAIN)
    predictions_train = model.predict_proba(X_TRAIN)
    AUPRC_train = average_precision_score(Y_TRAIN,predictions_train[:,1])
    range_train_predictions = max(predictions_train[:,1])-min(predictions_train[:,1])
    results_train[name]=(AUPRC_train_cv,AUPRC_train,range_train_predictions)
    predictions_test = model.predict_proba(X_TEST)
    AUPRC_test = average_precision_score(Y_TEST,predictions_test[:,1])
    results_test[name]=AUPRC_test,predictions_test
    
names=[]
values=[]
for key in results_test:
    names.append(key)
    values.append(results_test[key][0])

max_result = max(values)
for algorithm in list(results_test.keys()):
    if results_test[algorithm][0] < 0.7*max_result or ((results_train[algorithm][0]-results_test[algorithm][0])**2)**0.5 > 0.2 or results_train[algorithm][2]<0.8:
        with open(f'{out_dir}/deleted_algorithms.txt','a') as f:
            f.write(f'{algorithm}: test:{results_test[algorithm][0]}, train:{results_train[algorithm][0]},{results_train[algorithm][1]} range:{results_train[algorithm][2]}\n')
        results_test.pop(algorithm)
    else:
        with open(f'{out_dir}/remaining_algorithms.txt','a') as rem_f:
            rem_f.write(f'{algorithm}: test:{results_test[algorithm][0]}, train:{results_train[algorithm][0]},{results_train[algorithm][1]} range:{results_train[algorithm][2]}\n')

all_classifiers_preds_test=[]
all_classifiers_scores_test=[]
for key in results_test:
    all_classifiers_preds_test.append(results_test[key][1][:,1])
    all_classifiers_scores_test.append(results_test[key][0])

w_average=np.average(all_classifiers_preds_test,axis=0,weights=all_classifiers_scores_test)
w_av_AUPRC=average_precision_score(Y_TEST,w_average)
#%%
#feature importance
rdf=RandomForestClassifier(n_estimators=1000)
rdf.fit(X,Y)
importance=rdf.feature_importances_
sorted_feature_importance = sorted(zip(list(dataset.columns)[1:-1], importance), key=lambda x: x[1], reverse=True)
with open(f'{out_dir}/feature_importance.csv','w',newline='') as f:
    writer=csv.writer(f)
    writer.writerow(['feature','importance'])
    for feature in sorted_feature_importance:
        writer.writerow(feature)

predictions={}
algorithms_names = []
for name,model,X_TRAIN,X_TEST,Y_TRAIN,Y_TEST,all_X,all_data in updated_models:
    if name in results_test:
        model.fit(all_X,Y)
        model_proba=model.predict_proba(all_data)
        if max(model_proba[:,1])-min(model_proba[:,1]) >= 0.8:
            predictions[name]={}
            for i in range(len(model_proba)):
                predictions[name][array[i][0]]=model_proba[i][1]
            algorithms_names.append(name)
       
with open(f'{out_dir}/predictions.csv','w',newline='') as file:
    f_writer = csv.writer(file)
    header=['locus']
    weights=[]
    for name in algorithms_names:
        header.append(f'{name} (AUPRC {"%.3f" % results_test[name][0]})')
        weights.append(results_test[name][0])
    header.append(f'wVote (AUPRC {"%.3f" % w_av_AUPRC})')
    header.append('is_effector')
    f_writer.writerow(header)
    for locus in sorted(predictions['RDF'],key=predictions['RDF'].get,reverse=True):
        l=[locus]
        for name in algorithms_names:
            l.append(predictions[name][locus])
        l.append(np.average(l[1:len(algorithms_names)+1],weights=weights))
        if locus in list(yes.locus):
            l.append('yes')
        elif locus in list(no.locus):
            l.append('no')
        else:
            l.append('?')
        f_writer.writerow(l)

preds_df = pd.read_csv(f'{out_dir}/predictions.csv')

classifiers_scores = {f'{name} (AUPRC {"%.3f" % results_test[name][0]})':results_test[name][0] for name in results_test}
classifiers_scores[f'wVote (AUPRC {"%.3f" % w_av_AUPRC})'] = w_av_AUPRC
    
best_classifier = max(classifiers_scores,key= lambda key: classifiers_scores[key])
#if classifiers_scores[best_classifier]-w_av_AUPRC > 0.02: # if another classifier is significantly better than the vote, take it.
#    concensus_df = preds_df[['locus',f'{best_classifier}','is_effector']]
#    concensus_df.rename(columns = {best_classifier:f'likelihood to be effector by {best_classifier}'},inplace=True)
#    sorted_concensus = concensus_df.sort_values(by=f'likelihood to be effector by {best_classifier}',ascending=False)
#    sorted_concensus.to_csv(f'{out_dir}/concensus_predictions.csv',index=False)
#else: # vote is the best classifier or very close to it, so take it
concensus_df = preds_df[['locus',f'wVote (AUPRC {"%.3f" % w_av_AUPRC})','is_effector']]
sorted_concensus = concensus_df.sort_values(by=f'wVote (AUPRC {"%.3f" % w_av_AUPRC})',ascending=False)
sorted_concensus.rename(columns = {f'wVote (AUPRC {"%.3f" % w_av_AUPRC})':f'score (AUPRC on test set: {"%.3f" % w_av_AUPRC})'},inplace=True)
sorted_concensus.to_csv(f'{out_dir}/concensus_predictions.csv',index=False)

# figures

importance_f = f'{out_dir}/feature_importance.csv'
best_features =[]
with open(importance_f,'r') as in_f:
    reader = csv.reader(in_f)
    next(reader)
    for i in range(10):
        feature = next(reader)[0]
        best_features.append(feature)
f = features_file
#f = r'C:\Users\TalPNB2\Naama\Naama\effectors\Citrobacter_rodentium\revision\Citrobacter_features_with_addition.csv'
data = pd.read_csv(f)
labeled=data[data.is_effector!='?']
f_names = labeled.columns[1:-1]
labeled[f_names] = labeled[f_names].astype(float)
for i in range(len(best_features)):
    fig = plt.figure(figsize=(8,6))
    ax = sns.violinplot(x="is_effector", y=best_features[i], data=labeled)
    y_label = ' '.join(best_features[i].split('_'))
    if len(y_label)>45:
        y_lable_l = best_features[i].split('_')
        y_label = ' '.join(y_lable_l[:len(y_lable_l)//2])+'\n'+' '.join(y_lable_l[len(y_lable_l)//2:])
    ax.set(ylabel=y_label)
    fig.savefig(f'{out_dir}/{str(i)}.png')


df=pd.read_csv(importance_f)
fig=plt.figure(figsize=(10,5))
ax=sns.barplot(x='importance',y='feature',data=df.head(n=10))
ax.set_xlabel('')
ax.set_ylabel('',fontsize=16)
#ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
#fig.suptitle('Feature importance\n',y=1)
fig.subplots_adjust(left=0.5)
fig.tight_layout()

fig.savefig(f'{out_dir}/feature_importance.png')