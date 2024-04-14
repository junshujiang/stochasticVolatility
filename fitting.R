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
  ),
  list( 
    model = "rs",
    Y="stochvol.nig"
  ),
  list( 
    model = "rs",
    Y="stochvol.t"
  ),
  list( 
    model = "rs",
    Y="stochvolln"
  )
)
start_values=list(
  list( #1
    prior.std.dev.nominal=2,
    rspde.order=3,
    prior.range.nominal=1
    ),
  list(# 2
    prior.std.dev.nominal=2,
    rspde.order=5,
    prior.range.nominal=1
    ),
  list( #3
    prior.std.dev.nominal=2,
    rspde.order=7,
    prior.range.nominal=1
    ),
  list( #4
    prior.std.dev.nominal=20,
    rspde.order=3,
    prior.range.nominal=10
    ),
  list( #5
    prior.std.dev.nominal=20,
    rspde.order=5,
    prior.range.nominal=10
    ),
  list( #6
    prior.std.dev.nominal=20,
    rspde.order=7,
    prior.range.nominal=10
    ),
  list( #7
    prior.std.dev.nominal=50,
    rspde.order=3,
    prior.range.nominal=30
    ),
  list( #8
    prior.std.dev.nominal=50,
    rspde.order=5,
    prior.range.nominal=30
    ),
  list( #9
    prior.std.dev.nominal=50,
    rspde.order=7,
    prior.range.nominal=30
    )
)
configs_nums=length(families)*length(start_values)
configs=list()
config_i=1

for (family in families){
  for (start_value in start_values){
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
for (config_i in 1:configs_nums){
  config<-configs[[config_i]]
  family<-config$family
  startPoints<-config$start_value
  fit<-buildModel(date_timeindex,logreturn,startPoints,family,fit_times =10) 
  config$fit<-fit
  saveRDS(config, file.path(model_save_path, glue("fit_{config_i}.rds")))
  futile.logger::flog.info(glue("fit_{config_i} is saved"),name = "xx")
}

