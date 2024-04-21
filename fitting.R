source('helper.R')
library(readr)
library(glue)
source("helper.R")

max_loger_id<-get_max("Log")
expName<-"fitting"

exp_path<-glue("Log/{max_loger_id}.{expName}")
logger=get_logger(exp_path)

data_ori <- as.data.frame(read_csv("tmp/1180_vol.csv"))
logpt<-log(data_ori$Close)#
logreturn <-logpt[2:length(logpt)]-logpt[1:(length(logpt)-1)]
date <- data_ori$Date[1:(length(data_ori$Date)-1)]
date_timeindex<-as.numeric(as.POSIXct(date))
date_timeindex<-(date_timeindex - min(date_timeindex))/86400

num_start_value=30

families <- list(
  list( 
    model = "rough",
    Y="stochvol"
  ),
  list( 
    model = "rough",
    Y="stochvol.nig"
  ),
  list( 
    model = "rough",
    Y="stochvol.t"
  ),
  list( 
    model = "rough",
    Y="stochvolln"
  ),
  list( 
    model = "rs",
    Y="stochvol"
  )
)

configs_nums=length(families)*num_start_value
configs=list()
config_i=1

for (family in families){
  for (start_value_i in 1:num_start_value){

    order.start =sample(2:6, 1)
    range.start =runif(1,min=1,max=300)
    if (family$model=="rough" & family$Y=="stochvolln" ){
      std.start = rlnorm(1,-0.49184250617290765, 0.42746915578970635)  

      }

    else if (family$model=="rough" & family$Y=="stochvol.nig" ) {
      std.start = rlnorm(1,-0.28759736750416326, 0.36233095420105155) 
    }
    else if (family$model=="rough" & family$Y=="stochvol.t" ) {
      std.start = rlnorm(1,-0.3337091281534198, 0.3764506474038931) 
    }
    else if (family$model=="rough" & family$Y=="stochvol" ) {
      std.start = rlnorm(1,-0.426690, 0.4059762) 
    }
    start_value=list(rspde.order=order.start,prior.range.nominal=range.start,prior.std.dev.nominal=std.start) 
    configs[[config_i]]<-list(
      family=family,
      start_value=start_value
    )
    config_i<-config_i+1
  }
}
futile.logger::flog.info(glue("Total configures {configs_nums}"), name = "xx")
model_save_path=file.path(exp_path, 'fits')
if (!dir.exists(model_save_path)) {
    dir.create(model_save_path, recursive = TRUE)
  }
for (config_i in 1:length(configs)){

  config<-configs[[config_i]]
  family<-config$family
  startPoints<-config$start_value
  #fit<-buildModel(date_timeindex,logreturn,startPoints,family,fit_times =10)
  result <- try(buildModel(date_timeindex,logreturn,startPoints,family,fit_times =1) ,silent = TRUE)
  
  model=family$model
  family=family$Y
  rspde.order=startPoints$rspde.order
  prior.range.nominal=startPoints$prior.range.nominal
  prior.std.dev.nominal=startPoints$prior.std.dev.nominal

  futile.logger::flog.info(glue("model:{model},Y:{family}"),name = "xx")
  futile.logger::flog.info(glue("rspde.order:{rspde.order},prior.range.nominal:{prior.range.nominal},prior.std.dev.nominal:{prior.std.dev.nominal}"),name = "xx")
  if (inherits(result, "try-error")) {
	  futile.logger::flog.info("Error",name = "xx")
    } else {
        
  config$fit<-result
  saveRDS(config, file.path(model_save_path, glue("fit_{config_i}.rds")))
  futile.logger::flog.info(glue("fit_{config_i} is saved"),name = "xx")
}
}

