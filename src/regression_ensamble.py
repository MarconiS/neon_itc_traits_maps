#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 14:50:56 2020

@author: sergiomarconi
"""


import numpy as np
import pandas as pd

import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.ensemble import StackingClassifier
from sklearn.ensemble import BaggingClassifier
from mlxtend.classifier import StackingCVClassifier
from sklearn.linear_model import LogisticRegressionCV
from sklearn.experimental import enable_hist_gradient_boosting 
from sklearn.ensemble import HistGradientBoostingClassifier
from mlxtend.classifier import SoftmaxRegression
from sklearn.ensemble import GradientBoostingClassifier


from imblearn.ensemble import BalancedRandomForestClassifier
from imblearn.ensemble import RUSBoostClassifier

from sklearn.svm import SVC

data = pd.
#define models
rf = RandomForestRegression(random_state=0, oob_score = True, n_jobs = 2, 
                            n_estimators = 500, max_features = 'sqrt', criterion = 'entropy')
#clf3 = GaussianNB()

#gb = GradientBoostingClassifier( random_state=0)
#mlp = MLPClassifier(solver='lbfgs',random_state=0)
gb = HistGradientBoostingRegressor(random_state=0, max_iter = 1000, learning_rate = 0.1, 
                max_depth = 25, l2_regularization = 0.5)
bsvc = BaggingRegressor(base_estimator=SVR(), n_jobs = 2, 
                         oob_score = True, random_state=0)


# Initializing models
clf_bl = StackingCVRegressor(classifiers = [make_pipeline(StandardScaler(),rf), 
                                             make_pipeline(StandardScaler(),gb), 
                                             make_pipeline(StandardScaler(),clf)],
                          use_probas=True,
                          average_probas=False,
                          meta_classifier= LogisticRegressionCV())

params = {
 'meta_classifier__Cs': [0.1, 5, 10], 
 'meta_classifier__max_iter': [10000],
 }

grid = GridSearchCV(estimator=clf_bl, 
                    param_grid=params, 
                    cv=3,
                    refit=True)
grid.fit(X_res, y_res.taxonID.ravel())


clf_bl.fit(X_res, y_res.taxonID.ravel())
print(clf_bl.score(X_test, y_test['taxonID'].ravel()))
#rf_check = brf.fit(X_res, y_res.taxonID)

# #hipertune imbalanced models
# params = {'kneighborsclassifier__n_neighbors': [1, 5],
#           'randomforestclassifier__n_estimators': [10, 50],
#           'meta_classifier__C': [0.1, 10.0]}

# grid = GridSearchCV(estimator=sclf, 
#                     param_grid=params, 
#                     cv=3,
#                     refit=True)
# grid.fit(X, y)
# cv_keys = ('mean_test_score', 'std_test_score', 'params')

predict_an = clf_bl.predict_proba(X_test)
predict_an = pd.DataFrame(predict_an)



from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import f1_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
#from class_hierarchy import *
# format outputs for final evaluation
rf = rf.fit(X_res, y_res.taxonID)
taxa_classes = rf.classes_
colnames = np.append(['individualID', 'taxonID'], taxa_classes) #.tolist()
y_test.reset_index(drop=True, inplace=True)
predict_an.reset_index(drop=True, inplace=True)
eval_an = pd.concat([y_test, predict_an], axis=1)
eval_an.columns = colnames
#aggregate probabilities of each pixel to return a crown based evaluation. Using majority vote
eval_an = eval_an.groupby(['individualID', 'taxonID'], as_index=False).mean()

y_itc = eval_an['taxonID']
pi_itc = eval_an.drop(columns=['individualID', 'taxonID'])

