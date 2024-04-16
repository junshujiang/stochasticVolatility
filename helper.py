
import logging
import os
import numpy as np 



def get_max(dirpath)->int:
    """
    Author:
    Description:
        Incremental ID generator, used to generate incremental ID for expriment
    Args:
        dirpath:str, Log path
    Returns: 
        max+1:int, Incremental ID
    Raises:
    Example:
    Update:
    """
    array=np.array([int(x.split(".")[0]) for x in os.listdir(dirpath) if not x.startswith(".")])
    if len(array)==0:
        return 1
    else:
        return array.max()+1
    

dict_logger=dict()
def get_logger(log_path,log_name="xx",debug=False):
    """
    Author:
    Description:
        get logger object, so you can write log easily
    Args:
        log_path:str, Log path
        log_name: log_name: help you to identify your log file
        debug: bool, whether to print log in console
    Returns: 
        logger
    Raises:
    Example:
    Update:
    """
    global dict_logger
    
    if log_path+log_name not in dict_logger:
        logger = logging.getLogger(log_path+log_name)
        dict_logger[log_path+log_name]=logger
        logger.setLevel(logging.INFO)  


        if log_path=="":
            new_log_path = os.path.abspath(".") + '/Logs/'
        else:
            new_log_path=log_path+"/"
        if not os.path.exists(new_log_path):
            os.mkdir(new_log_path)

        logfile = new_log_path+ log_name+'.log'
        fh = logging.FileHandler(logfile, mode='a')
        fh.setLevel(logging.DEBUG)  
        formatter = logging.Formatter("%(asctime)s- %(levelname)s: %(message)s")
        fh.setFormatter(formatter)


        if debug:
            ch=logging.StreamHandler()
            ch.setLevel(logging.INFO)
            ch.setFormatter(formatter)
            logger.addHandler(ch)

        logger.addHandler(fh)
        
    return dict_logger[log_path+log_name]

