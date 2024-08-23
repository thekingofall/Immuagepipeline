from importlib.resources import open_binary
import xgboost as xgb
import numpy as np
import scipy
import pandas as pd


import joblib

from xgboost import plot_importance
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split #nolint
from numpy import absolute

from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedKFold
import os


from optparse import OptionParser
import time,datetime

parser = OptionParser(usage="usage: %prog -f readFolde [-p <int>] [-e <int>] ", version="%prog 1")
parser.add_option("-f","--fo",
        action="store", 
        dest="fo",
        help="the file about gene expression matrix.default:featureold.csv")
parser.add_option("-c","--coldata",
        action="store", 
        dest="coldata",
        help="the file about metadata:default is 'coldata.csv'")
parser.add_option("--newfeature",
        action="store_true",    
        dest="newfeature",
        help="the file about gene expression matrix.default:featurenew.csv")
parser.add_option("--newcol",
        action="store",
        dest="newcol",
        help="the file about gene expression matrix.default:newcol.csv")



parser.add_option("-m","--modelchoose",
        action="store", 
        dest="modelchoose",
        default="xgb",
        help="the model file to save.default is 'xgb'")

parser.add_option("-l","--modelload",
        action="store", 
        dest="modelload",
        help="the model file to load.default is 'xgb'")

parser.add_option("-p","--project",
        action="store",
        dest="project",
        default="PBMC",
        help="the project name.default is 'PBMC'")



parser.add_option("-w","--ways",
        action="store", 
        dest="ways",
        default="train",
        help="the way of the method is 'predict' or 'trainpre''train'.default is 'train'")
parser.add_option("-e","--eval",
        action="store", 
        dest="eval",
        help="the file to evaluate.default is 'eval.csv'")

parser.add_option("--rootdir",
        action="store",
        dest="rootdir",
        help="the root directory of the project.default is '.'")


(options, args) = parser.parse_args()


def train_model(X, y, params, num_boost_round=100, early_stopping_rounds=10,mode="xgb"):
    if mode == "xgb":
    # train model
        xgtrain = xgb.DMatrix(X, label=y)
        model = xgb.XGBRegressor(params, xgtrain, num_boost_round, early_stopping_rounds=early_stopping_rounds)
        return model
    elif mode == "auto":
        from autogluon.tabular import TabularPredictor

        df_train = X
        df_train["age"] = y.values
        predictor= TabularPredictor(label ='age').fit(train_data = df_train, verbosity = 2,presets='best_quality')
        return predictor



if __name__ == '__main__':
    (options, args) = parser.parse_args()

    if options.ways == 'train':

        if not options.fo  or not options.coldata :
            parser.print_help()
            exit(0)
        else:
            featureold = options.fo
            coldata = options.coldata
            print(featureold)
            print(coldata)
            featureolds=pd.read_csv(featureold,sep='\t',header=0)
            coldatas=pd.read_csv(coldata,sep='\t',header=0)
            print(featureolds.shape)
            print(coldatas.head())
            
            coldatas.index=coldatas['name']
            features=featureolds.T
            print(coldata)
            print(features.shape)
            target=coldatas["age"]
            if list(features.index) != list(coldatas.index):
                print('the index of feature and coldata are not equal')
                exit(0)
            else:
                print('the index of feature and coldata are equal')
                x_train,x_test,y_train,y_test=train_test_split(features,target,test_size=0.2,shuffle=True)
                print(x_train.shape)
                #parms={"learning_rate":0.1,"n_estimators":100,"max_depth":10,"min_child_weight":1,"gamma":0,"subsample":0.8,"colsample_bytree":0.8,"scale_pos_weight":1,"objective":"reg:linear","eval_metric":"rmse"}
               # model=train_model(x_train,y_train,parms,num_boost_round=100,early_stopping_rounds=10,mode="xgb")
                model = xgb.XGBRegressor(n_estimators=1000, max_depth=10, eta=0.1, subsample=0.9, colsample_bytree=0.8)
                model.fit(x_train,y_train)
                y_pred=model.predict(x_test)
                print(y_pred)
                print(y_test)
                modelsavedir=os.path.join(os.getcwd(),'A06_model')
                os.system("mkdir -p A06_model")
                endtime =str(datetime.datetime.now()).split(' ')[0]
                plt.figure(figsize=(8, 6))
                plot_importance(model,max_num_features=10)
            
                plt.savefig(modelsavedir+"/"+endtime+'featureimportance.png', dpi=500)
                
                
                if options.modelchoose=="xgb":

                    joblib.dump(model,modelsavedir +"/xgb.dat")
                    print("xgb model saved")

                elif options.modelchoose=="auto":
                    model.save(modelsavedir +'/auto.model')
                    print("auto model saved")

    elif options.ways == 'predict':
        if not options.modelload and not options.modelchoose:
            parser.error("trained  model must be given and what kinds of module must be given:-w ,-m,-l")
            exit(0)
        else:
            featurenew = options.fo
            loaded_model = joblib.load(options.modelload)
            if options.coldata :
                coldata = pd.read_csv(options.coldata,sep='\t',header=0)
                print(list(coldata["age"]))
            featureolds=pd.read_csv(featurenew,sep='\t',header=0)
            features=featureolds.T
            print(features.shape)
            #model = xgb.XGBRegressor(n_estimators=1000, max_depth=10, eta=0.1, subsample=0.9, colsample_bytree=0.8)
        
            y_pred=loaded_model.predict(features)
            print(y_pred)
            endtime =str(datetime.datetime.now()).split(':')[0].split(' ')[0]
            print(endtime)
            writename=options.project+"_"+endtime+str(options.modelchoose)+".csv"
            with open(writename,'w') as f:
                f.write("name,age\n")
                for i in range(len(y_pred)):
                    f.write(str(features.index[i])+","+str(y_pred[i])+"\n")

    elif options.ways =='train_pre':
        if not options.fo  or not options.coldata :
            parser.print_help()
            exit(0)
        else:
            featureold_arg = options.fo
            coldata_arg = options.coldata
            newcoldata_arg=options.newcol
            print(newcoldata_arg)
     
            features=pd.read_csv(featureold_arg,sep='\t',header=0)
            coldatas=pd.read_csv(coldata_arg,sep='\t',header=0)
            newcoldatas=pd.read_csv(newcoldata_arg,sep='\t',header=0)
         
            print(coldatas.head())
            
            coldatas.index=coldatas['name']
            print(coldatas.index)
            featureolds=features.T
            print(coldatas)
            print(featureolds.head())
            print(featureolds.shape)
            print(features.shape)
            target=coldatas["age"]
            print(coldatas)
            print(featureolds.index)
            if list(featureolds.index) != list(coldatas.index):
                print('the index of feature and coldata are not equal')
                exit(0)
            else:
                print('the index of feature and coldata are equal')
                newname=newcoldatas["name"]
                print(newcoldatas)
                feature_new=featureolds.loc[newcoldatas["name"],:]
                featureolds.drop(newcoldatas["name"],axis=0,inplace=True)
                feature_oldfortrain=featureolds
                print(feature_oldfortrain.head)
                print(feature_oldfortrain.shape)
                print(feature_new.shape)
                coldatas.drop(newcoldatas["name"],axis=0,inplace=True)
                target_new=coldatas["age"]
                print(target_new.shape)

                x_train,x_test,y_train,y_test=train_test_split(feature_oldfortrain,target_new,test_size=0.2,shuffle=True)
                print(x_train.shape)
                #parms={"learning_rate":0.1,"n_estimators":100,"max_depth":10,"min_child_weight":1,"gamma":0,"subsample":0.8,"colsample_bytree":0.8,"scale_pos_weight":1,"objective":"reg:linear","eval_metric":"rmse"}
               # model=train_model(x_train,y_train,parms,num_boost_round=100,early_stopping_rounds=10,mode="xgb")
                model = xgb.XGBRegressor(n_estimators=1000, max_depth=10, eta=0.1, subsample=0.9, colsample_bytree=0.8)
                model.fit(x_train,y_train)
                # print(model.feature_importances_)
                # print(model.best_score_)
                # print(model.best_iteration_)
                # print(model.best_ntree_limit_)

                y_pred=model.predict(x_test)
                plot_importance(model,max_num_features=10)
                plt.savefig(os.getcwd()+"/featureimportance.png", dpi=500)

                print(y_pred)
                print(y_test)
                # from sklearn.metrics import mean_squared_error
                # mse_score = mean_squared_error(y_test, y_pred)
                # print("Mse score is:", mse_score)
                # print(model.score(x_test,y_test))
                from sklearn.metrics import r2_score #直接调用库函数进行输出R2
                from sklearn.metrics import mean_absolute_error
                from sklearn.metrics import mean_squared_error


                R2=r2_score(y_test.values,y_pred)
                MAE=mean_absolute_error(y_test.values,y_pred)
                MSE=mean_squared_error(y_test.values,y_pred)

                print("R2 is",r2_score(y_test.values,y_pred))
                print("MAE is ",mean_absolute_error(y_test.values,y_pred))
                print("MSE is ",mean_squared_error(y_test.values,y_pred))
                data_metrics = pd.DataFrame([{"R2":R2,"MAE":MAE,"MSE":MSE}])
                data_metrics.to_csv(options.project+"_"+str(options.modelchoose)+"_Mteric.csv",index=False)  


                new_pred=model.predict(feature_new)
                print("新增人员信息","\n",newcoldatas)
                print("\n")

                print("测试集预测结果")
                data_test=pd.DataFrame({"name":y_test.index," real_age":y_test.values,"predict_age":y_pred})
                print(data_test)
                # print("测试集预测年龄：",y_pred)
                # print("测试集实际年龄：",y_test)

                print("\n")
                print("新增预测年龄：",new_pred)
                print("新增实际年龄：",newcoldatas["age"].values)
                

                
            

                os.system("rm -rf "+options.project+"_"+str(options.modelchoose)+"_testpredict_results.csv")
                results=options.project+"_"+str(options.modelchoose)+"_testpredict_results.csv"

                data_test.to_csv(results,index=False,header=True,sep=',')
                # results.write(data_test+"\n") 
                
              
                # results.write(y_pred,y_test+"\n")
               
          
                # results.write(newname,new_pred,newcoldatas["age"])

                print("测试集预测年龄：",y_pred)
                print("测试集实际年龄：",y_test)
                print("新增预测年龄：",new_pred)
                print("新增实际年龄：",newcoldatas["age"])
                modelsavedir=os.path.join(options.rootdir,'A06_model')
                os.system("mkdir -p "+modelsavedir)
                if options.modelchoose=="xgb":
                    joblib.dump(model,modelsavedir +"/xgb.dat")
                    print("xgb model saved")



        









