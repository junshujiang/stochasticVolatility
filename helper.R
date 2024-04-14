library(INLA) ##
library(rSPDE)
library(readr)
library(fmesher)
library(inlabru)
library(futile.logger)


buildModel <- function(date_timeindex,logreturn,startPoints,family,fit_times=1) {
    # Note: data_timeindex should include all the points in the time series, not just the training set
       
    # Author:
    # Description:

    # Args: 
    #     date_timeindex: build mesh
    #     logreturn: time series
    #     startPoints: startPoints for params 
    #     family: 
    # Returns:
    # Raises: 
    # Example:
    # Update: 
    result_total=list()
    mesh <- fm_mesh_1d(date_timeindex) 
    data_length=length(logreturn)
    data = list(logreturn = logreturn, times = date_timeindex[1:data_length])
    if (family$model == "rough"){

    rspde_model <- rspde.matern(mesh = mesh,prior.std.dev.nominal=startPoints$prior.std.dev.nominal,
    rspde.order=startPoints$rspde.order, prior.range.nominal=startPoints$prior.range.nominal)
    formula <- logreturn ~ f(times, model = rspde_model)



    #result_str<-summary(result_fit)
    }
    else{
    rspde_model1 <- rspde.matern(mesh = mesh,prior.std.dev.nominal=startPoints$prior.std.dev.nominal,
    rspde.order=startPoints$rspde.order, prior.range.nominal=startPoints$prior.range.nominal)
    
    rspde_model2 <- rspde.matern(mesh = mesh,prior.std.dev.nominal=startPoints$prior.std.dev.nominal,
    rspde.order=startPoints$rspde.order, prior.range.nominal=startPoints$prior.range.nominal)
    formula <- logreturn ~ Intercept(1) + rough(times, model = rspde_model1) + smooth(times, model = rspde_model2)
    }


    fit_nums=0
    while(fit_nums<fit_times){
        result=list()
        fit <- bru(formula, family=family$Y, data = data)
        if (family$model == "rough"){
            result_fit<- rspde.result(fit, "f", rspde_model, parameterization="matern")
            result$fit<-fit
            result$spde_fit<-result_fit

        }
        else{
        result_fit1<- rspde.result(fit, "rough", rspde_model1, parameterization="matern")
        result_fit2<- rspde.result(fit, "smooth", rspde_model2, parameterization="matern")
        result$fit<-fit
        result$rough_fit<-result_fit1
        result$smooth_fit<-result_fit2

        }
        fit_nums=fit_nums+1
        result_total[[fit_nums]]<-result
    }


    
    return (result_total)
    }
    


parse_result<-function(fit){
    likelihoods=list()
    df <- data.frame()

    for( f_i in seq(length(fit))){

        f=fit[[f_i]]
        likelihoods=append(likelihoods,summary(f$fit)[2]$inla$mlik[1])
        if (family$model=="rough"){
            result_df=summary(f$spde_fit)

        }else{
            result_rough_df=summary(f$rough_fit)
            rownames(result_rough_df) <- paste0("rough_", rownames(result_rough_df))
            result_smooth_df=summary(f$smooth_fit)
            rownames(result_smooth_df) <- paste0("smooth_", rownames(result_smooth_df))
            result_df<-rbind(result_rough_df,result_smooth_df)
        }
        f_i_index<-f_i
        rownames(result_df) <- paste0(glue("e{f_i_index}"), rownames(result_df))
        df=rbind(df,result_df)

}
return(list(df=df,like=likelihoods))
}



fit_a_model_ar<-function(logreturn,date_timeindex, family="stochvol.nig"){
    data = list(logreturn = logreturn, times = date_timeindex)
    formula<-logreturn~f(times,model="ar1")
    fit=bru(formula,family=family,data=data)
    return (fit)
    }

fit_a_model_ou<-function(logreturn,date_timeindex,nu=0.5,family="stochvol.nig"){

    mesh <- fm_mesh_1d(date_timeindex) 
    data_length=length(logreturn)
    data = list(logreturn = logreturn, times = date_timeindex[1:data_length])
    if( is.null(nu)){
    rspde_model <- rspde.matern(mesh = mesh)
    }
    else{
    rspde_model <- rspde.matern(mesh = mesh,nu=nu)
    }


    formula <- logreturn ~ f(times, model = rspde_model)
    fit <- bru(formula, family=family, data = data)
    return (fit)
    }

fit_a_model_ls<-function(logreturn,datepredict_ar_rolling_timeindex,family="stochvol.nig"){

    mesh <- fm_mesh_1d(date_timeindex)
    data_length=length(logreturn)
    data = list(logreturn = logreturn, times = date_timeindex[1:data_length])
    rspde_model1 <- rspde.matern(mesh = mesh)
    rspde_model2 <- rspde.matern(mesh = mesh, nu = 1)
    formula <- logreturn ~ Intercept(1) + rough(times, model = rspde_model1) + smooth(times, model = rspde_model2)
    fit <- bru(formula, family=family, data = data)
    return (list(fit=fit,rspde_model1=rspde_model1,rspde_model2=rspde_model2))
}

predict_ar<-function(train_logreturn, train_date_timeindex, test_logreturn, test_date_timeindex,windowSize=10){


    results=c()
    batch_num=ceiling(length(test_logreturn)/windowSize)
    all_logreturn <- c(train_logreturn,test_logreturn)
    all_date_timeindex <- c(train_date_timeindex,test_date_timeindex)

    for( i in 1:batch_num){
        fit<-fit_a_model_ou(train_logreturn,date_timeindex, family="stochvol",nu=0.5)
        test_end=min(length(train_logreturn)+(windowSize*(i)),length(train_date_timeindex))
        tmp_logreturn<-all_logreturn[1:test_end]
        tmp_date_timeindex<-all_date_timeindex[1:test_end]
        tmp_logreturn[(test_end-windowSize+1):test_end]=NA
        test_data = list(logreturn = tmp_logreturn, times = tmp_date_timeindex)
        field_pred <- predict(fit, test_data, ~exp(0.5*(Intercept + f)),n.samples=100)
        results=c(results,field_pred[(test_end-windowSize+1):test_end,1])


    }

    return (results)

}


## rolling predict future 10 and average 
predict_ar_rolling<-function(logreturn, timeindex, test_start,family="stochvol",windowSize=10){
    results=c()
    for( i in seq(test_start,length(logreturn)-windowSize)){
        train_logreturn<-logreturn[1:i]
        fit<-fit_a_model_ou(train_logreturn,timeindex, family=family,nu=0.5)
        tmp_date_timeindex<-timeindex[(i+1):(i+windowSize)]
        test_data = list(times = tmp_date_timeindex)
        field_pred <- predict(fit, test_data, ~exp(0.5*(Intercept + f)),n.samples=100)
        next_seq=mean(field_pred[,1])
        results=c(results,next_seq)

    }

    return (results)

}


## rolling predict future 10 and average 
predict_ls_rolling<-function(logreturn, timeindex, test_start,family="stochvol",windowSize=10){
    results=c()
    for( i in seq(test_start,length(logreturn)-windowSize)){
        train_logreturn<-logreturn[1:i]
        result<-fit_a_model_ls(train_logreturn,timeindex, family=family)
        fit=result$fit
        tmp_date_timeindex<-timeindex[(i+1):(i+windowSize)]
        test_data = list(times = tmp_date_timeindex)
        field_pred <- predict(fit, test_data, ~exp(0.5*(Intercept + rough+smooth)),n.samples=100)
        next_seq=mean(field_pred[,1])
        results=c(results,next_seq)

    }

    return (results)

}





get_logger <- function(log_path, log_name = "xx", debug = FALSE) {
  # 确保日志路径存在
  if (log_path == "") {
    log_path <- file.path(getwd(), 'Logs')
  } else {
    log_path <- file.path(log_path)
  }
  
  if (!dir.exists(log_path)) {
    dir.create(log_path, recursive = TRUE)
  }

  # 创建日志路径和文件名
  log_file <- file.path(log_path, paste0(log_name, '.log'))

  # 创建或获取日志器
  logger <- futile.logger::flog.logger(name = log_name)
  futile.logger::flog.threshold(futile.logger::INFO, name = log_name)
  futile.logger::flog.appender(appender.file(log_file), name = log_name)

  # 如果启用调试模式，则也在控制台输出日志
  if (debug) {
    futile.logger::flog.appender(appender.console(), name = log_name)
  }

  return(logger)
}
get_max <- function(dirpath) {
  # 检查目录中的文件
  filenames <- list.files(path = dirpath, full.names = FALSE)
  # 过滤出以数字开始且后接点的文件名
  numeric_files <- grep("^[0-9]+\\.", filenames, value = TRUE)
  # 从文件名中提取前面的整数部分
  ids <- as.integer(sub("\\..*$", "", numeric_files))  # 删除第一个点及其后的所有内容，并转换为整数
  # 计算最大 ID 并增加 1，如果没有文件则返回 1
  if (length(ids) == 0) {
    return(1)
  } else {
    return(max(ids) + 1)
  }
}


# # 首先，创建或获取日志器
# logger <- get_logger("Log/2.name", log_name = "example_log", debug = TRUE)

# # 使用日志器记录信息
# futile.logger::flog.info("This is an informational message", name = "example_log")

# # 记录警告
# futile.logger::flog.warn("This is a warning message", name = "example_log")

# # 记录错误
# futile.logger::flog.error("This is an error message", name = "example_log")

# # 记录严重错误
# futile.logger::flog.fatal("This is a fatal error message", name = "example_log")