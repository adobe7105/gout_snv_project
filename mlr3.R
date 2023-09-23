#特征工程snv############################################################################
library("mlr3")
library("mlr3verse")
library("mlr3extralearners")
library("data.table")
library("ggplot2")
library("caret")
set.seed(1567)
lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")

set.seed(1567)
snv_trainIndex <- createDataPartition(snv188t$target_snv, p = .8, list = FALSE, times = 1)
snv_find2_find=snv188t[snv_trainIndex,]
snv_find2_control=snv188t[-snv_trainIndex,]

a2.lrn_ranger_importance1   = lrn("classif.ranger",id="Ranger", importance = "impurity",predict_type="prob",num.threads=32)

snv_gout_0w <- as_task_classif(snv140t,target="target_snv")
set.seed(1567)
instance_ranger = fsi(
  task = snv_gout_0w ,
  learner = a2.lrn_ranger_importance1  ,
  resampling = rsmp("cv", folds = 10),
  measure = msr("classif.ce"),
  #terminator = trm("evals", n_evals = 50),
  terminator = trm("none")
)

set.seed(1567)
fselector = fs("rfe",recursive=T,subset_sizes=c(20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2))
fselector$optimize(instance_ranger)


library(viridisLite)
library(mlr3misc)
data_plot= as.data.table(instance_ranger$archive)
data_plot[,n:=map_int(importance,length)]
data_plot=filter(data_plot,n<=20)
ggplot(data_plot, aes(x =n, y = classif.ce)) +
  geom_line(
    color = viridis(1, begin = 0.5),
    linewidth = 1) +
  geom_point(
    fill = viridis(1, begin = 0.5),
    shape = 21,
    size = 3,
    stroke = 0.5,
    alpha = 0.8) +  
  geom_vline(
    xintercept=5,
    #data_plot[classif.ce== min(classif.ce)]$n
    colour="red",
    linetype=3
  )+
  xlab("Number of Features randomforest") +
  scale_x_reverse() +
  theme_minimal()


p1=filter(data_plot,n==5)$features[[1]]
p1


snv_find=snv140t[,c("target_snv",p1)]
snv_control=snv48t[,c("target_snv",p1)]


snv_0w_gout_gout=as_task_classif(snv_find,target = "target_snv")
snv_0w_gout_control=as_task_classif(snv_control,target = "target_snv")

res=rsmp("cv",folds=10)

aa=resample(snv_0w_gout_control,a2.lrn_ranger_importance1,resampling = res)
aa$aggregate(msr("classif.auc"))
#模型调参############################################################
a2.lrn_ranger_importance1   = lrn("classif.ranger",id="Ranger", importance = "impurity",predict_type="prob",num.threads=16)
a2.lrn_xgboost_importance2  = lrn("classif.xgboost",id="Xgboost",predict_type="prob",nthread=64,nrounds=1000,early_stopping_rounds=500,
                                  early_stopping_set="test",eval_metric="error")
a2.lrn_lightgbm_importance3 = lrn("classif.lightgbm",id="Lightgbm",predict_type="prob",num_threads=16)

res=rsmp("cv",folds=10)

snv_find=snv140t[,c("target_snv",p1)]
snv_control=snv48t[,c("target_snv",p1)]

snv_0w_gout_gout=as_task_classif(snv_find,target = "target_snv")
snv_0w_gout_control=as_task_classif(snv_control,target = "target_snv")

snv_find_num=snv_find
snv_find_num[,2:6]=lapply(snv_find_num[,2:6],as.numeric)
snv_find_num[,2:6]=snv_find_num[,2:6]-1
snv_control_num=snv_control
snv_control_num[,2:6]=lapply(snv_control_num[,2:6],as.numeric)
snv_control_num[,2:6]=snv_control_num[,2:6]-1

task_num_find=as_task_classif(snv_find_num,target = "target_snv")
task_num_control=as_task_classif(snv_control_num,target = "target_snv")

tune_space_xgboost=lts("classif.xgboost.default")
as.data.table(tune_space_xgboost)
##xgboost##############################
lrn_xgboost=lrn("classif.xgboost",id="Xgboost",predict_type="prob",nthread=64,
                nrounds=1000,early_stopping_rounds=500,
                early_stopping_set="test",eval_metric="error",
                eta=to_tune(1e-04,1,logscale=TRUE),
                max_depth =to_tune(1e+00,20,logscale=FALSE),
                colsample_bytree =to_tune(1e-01,1,logscale=FALSE),
                lambda=to_tune(1e-03,1000,logscale=TRUE),
                alpha=to_tune(1e-03,1000,logscale=TRUE) 
                )

instrance_xgboost=tune(
  method = tnr("random_search",batch_size=10),
  task=task_num_find,
  learner=lrn_xgboost,
  resampling = rsmp("cv",folds=10),
  measures = msr("classif.auc"),
  term_evals = 50,
  callbacks = clbk("mlr3tuning.early_stopping")
)

instrance_xgboost$result_learner_param_vals
learner_xgboost=lrn("classif.xgboost",predict_type="prob")
learner_xgboost$param_set$values=instrance_xgboost$result_learner_param_vals

ax=resample(task_num_find,learner_xgboost,resampling = res)
set.seed(3456)
ax$aggregate(msr("classif.auc"))
autoplot(ax,type="roc")

#svm##############################

lrn_svm=lrn("classif.svm",
            type="C-classification",
            kernel="radial",
            cost=to_tune(1e-1,1e5),
            gamma=to_tune(1e-1,1),
            predict_type="prob"
            )
instance_svm=ti(
  task = task_num_find,
  learner=lrn_svm,
  resampling = res,
  measures = msr("classif.auc"),
  terminator =trm("none")
)

tuner_svm=tnr("grid_search",resolution=10,batch_size=10)
tuner_svm$optimize(instance_svm)

learner_svm=lrn("classif.svm",predict_type="prob")
learner_svm$param_set$values=instance_svm$result_learner_param_vals


as=resample(task_num_control,learner_svm,resampling = res)
set.seed(3456)
as$aggregate(msr("classif.auc"))
autoplot(as,type="roc")
#ranger###################################################################################

lrn_ranger=lrn("classif.ranger",
            alpha=to_tune(1e-1,1),
            max.depth =to_tune(1e+00,20,logscale=FALSE),
            mtry=to_tune(2,5),
            num.trees=to_tune(10,1000),
            predict_type="prob"
)
instance_ranger=ti(
  task = task_num_find,
  learner=lrn_ranger,
  resampling = res,
  measures = msr("classif.auc"),
  terminator =trm("none")
)

tuner_ranger=tnr("grid_search",resolution=5,batch_size=10)
tuner_ranger$optimize(instance_ranger)

learner_ranger=lrn("classif.ranger",predict_type="prob")
learner_ranger$param_set$values=instance_ranger$result_learner_param_vals


ar=resample(task_num_find,learner_ranger,resampling = res)
set.seed(3456)
ar$aggregate(msr("classif.auc"))
autoplot(ar,type="roc")
#######################################################################
design=benchmark_grid(tasks = task_num_control,
                      learners=lrns(c("classif.xgboost","classif.svm"),
                                    predict_type="prob"),
                      resamplings = res
                        )
set.seed(3456)
bmr=benchmark(design)
autoplot(bmr,type="roc")
design2=benchmark_grid(tasks = task_num_gout,
                       learners=lrns(c(learner_ranger,learner_xgboost,learner_svm)),
                       predict_type="prob",
                       resamplings = res
)
set.seed(3456)
bmr2=benchmark(design2)
autoplot(bmr2,type="roc")

#特征数量#######################################################################
library("mlr3")
library("mlr3verse")
library("mlr3extralearners")
library("data.table")
library("ggplot2")
library("caret")
set.seed(3456)
lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")

data_trainIndex <- createDataPartition(data_188$target, p = .8,  list = FALSE,times = 1)
head(data_trainIndex)
data_find2_find=data_188[data_trainIndex,]
data_find2_control=data_188[-data_trainIndex,]

a2.lrn_ranger_importance1   = lrn("classif.ranger",id="Ranger", importance = "impurity",predict_type="prob",num.threads=16)

data_gout_0w <- as_task_classif(data_188,target="target")
set.seed(3456)
instance_ranger2 = fsi(
  task = data_gout_0w ,
  learner = a2.lrn_ranger_importance1  ,
  resampling = rsmp("cv", folds = 10),
  measure = msr("classif.ce"),
  #terminator = trm("evals", n_evals = 50),
  terminator = trm("none")
)

set.seed(3456)
fselector = fs("rfe",recursive=T,subset_sizes=c(20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2))
fselector$optimize(instance_ranger2)


library(viridisLite)
library(mlr3misc)
data_plot= as.data.table(instance_ranger2$archive)
data_plot[,n:=map_int(importance,length)]
data_plot=filter(data_plot,n<=20)
ggplot(data_plot, aes(x =n, y = classif.ce)) +
  geom_line(
    color = viridis(1, begin = 0.5),
    linewidth = 1) +
  geom_point(
    fill = viridis(1, begin = 0.5),
    shape = 21,
    size = 3,
    stroke = 0.5,
    alpha = 0.8) +  
  geom_vline(
    xintercept=8,
    #data_plot[classif.ce== min(classif.ce)]$n
    colour="red",
    linetype=3
  )+
  xlab("Number of Features randomforest") +
  scale_x_reverse() +
  theme_minimal()


p2=filter(data_plot,n==8)$features[[1]]
p2


data_find=data_find[,c("target",p2)]
data_control=data_test[,c("target",p2)]

data_0w_gout_gout=as_task_classif(data_find,target = "target")
data_0w_gout_control=as_task_classif(data_control,target = "target")


res=rsmp("cv",folds=10)

aa=resample(data_0w_gout_gout,a2.lrn_ranger_importance1,resampling = res)
set.seed(3456)
aa$aggregate(msr("classif.auc"))
##XG#########################################################################
lrn_xgboost2=lrn("classif.xgboost",id="Xgboost",predict_type="prob",nthread=64,
                nrounds=1000,early_stopping_rounds=500,
                early_stopping_set="test",eval_metric="error",
                eta=to_tune(1e-04,1,logscale=TRUE),
                max_depth =to_tune(1e+00,20,logscale=FALSE),
                colsample_bytree =to_tune(1e-01,1,logscale=FALSE),
                lambda=to_tune(1e-03,1000,logscale=TRUE),
                alpha=to_tune(1e-03,1000,logscale=TRUE) 
)

instrance_xgboost2=tune(
  method = tnr("random_search",batch_size=10),
  task=data_0w_gout_gout,
  learner=lrn_xgboost2,
  resampling = rsmp("cv",folds=10),
  measures = msr("classif.auc"),
  term_evals = 50,
  callbacks = clbk("mlr3tuning.early_stopping")
)

instrance_xgboost2$result_learner_param_vals
learner_xgboost2=lrn("classif.xgboost",predict_type="prob")
learner_xgboost2$param_set$values=instrance_xgboost2$result_learner_param_vals

ax2=resample(data_0w_gout_control,learner_xgboost2,resampling = res)
set.seed(3456)
ax2$aggregate(msr("classif.auc"))
autoplot(ax2,type="roc")
##SVM###############################################################################

lrn_svm2=lrn("classif.svm",
            type="C-classification",
            kernel="radial",
            cost=to_tune(1e-1,1e5),
            gamma=to_tune(1e-1,1),
            predict_type="prob"
)
instance_svm2=ti(
  task = data_0w_gout_gout,
  learner=lrn_svm2,
  resampling = res,
  measures = msr("classif.auc"),
  terminator =trm("none")
)

tuner_svm2=tnr("grid_search",resolution=10,batch_size=10)
tuner_svm2$optimize(instance_svm2)

learner_svm2=lrn("classif.svm",predict_type="prob")
learner_svm2$param_set$values=instance_svm2$result_learner_param_vals


as2=resample(data_0w_gout_gout,learner_svm2,resampling = res)
set.seed(3456)
as2$aggregate(msr("classif.auc"))
autoplot(as2,type="roc")



##RF###############################################################################

lrn_ranger2=lrn("classif.ranger",
               alpha=to_tune(1e-1,1),
               max.depth =to_tune(1e+00,20,logscale=FALSE),
               mtry=to_tune(2,5),
               num.trees=to_tune(10,1000),
               predict_type="prob"
)
instance_ranger2=ti(
  task =data_0w_gout_gout,
  learner=lrn_ranger2,
  resampling = res,
  measures = msr("classif.auc"),
  terminator =trm("none")
)

tuner_ranger2=tnr("grid_search",resolution=5,batch_size=10)
tuner_ranger2$optimize(instance_ranger2)

learner_ranger2=lrn("classif.ranger",predict_type="prob")
learner_ranger2$param_set$values=instance_ranger2$result_learner_param_vals


ar2=resample(data_0w_gout_gout,learner_ranger2,resampling = res)
set.seed(3456)
ar2$aggregate(msr("classif.auc"))
autoplot(ar2,type="roc")


































