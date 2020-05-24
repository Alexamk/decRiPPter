# License: GNU Affero General Public License v3 or later

import logging
import os

class LoggerOperator():
    
    def __init__(self):
        # Have the logger handlers been set
        self.handlers_are_set = False
        
        # Which loggers has this operator returned?
        self.loggers_returned = []
        
    def return_logger(self, logname, include_thread):
        logger = logging.getLogger(logname)
        if self.handlers_are_set:
            self.setup_logger(logger, include_thread)
        self.loggers_returned.append((logger, include_thread))
        return logger
        
    def setup_logger(self, logger, include_thread):
        if include_thread:
            formatter = self.formatter_threaded
        else:
            formatter = self.formatter

        fh = logging.FileHandler(self.logfile)
        fh.setLevel(self.loglevel)
        fh.setFormatter(formatter)
        
        ch = logging.StreamHandler()
        ch.setLevel(self.loglevel)
        ch.setFormatter(formatter)
        
        logger.addHandler(fh)
        logger.addHandler(ch)
        
        logger.setLevel(self.loglevel)
    
    def update_loggers(self):
        for logger, include_thread in self.loggers_returned:
            self.setup_logger(logger, include_thread)
    
    def set_handlers(self, settings):
        self.loglevel = get_loglevel(settings)
        self.logfile = get_logfile(settings)
        self.formatter = logging.Formatter( '%(levelname)s - %(asctime)s - %(name)s - %(message)s')
        self.formatter_threaded = logging.Formatter( '%(levelname)s - %(asctime)s - %(name)s - %(threadName)s - %(message)s')
        self.handlers_are_set = True

    def startup(self, settings):
        self.set_handlers(settings)
        self.update_loggers()
        
def get_loglevel(settings):
    if settings['verbosity'] == 3:
        return logging.DEBUG
    if settings['verbosity'] == 2:
        return logging.INFO
    if settings['verbosity'] == 1:
        return logging.WARNING
    if settings['verbosity'] == 0:
        return logging.ERROR
    
def get_logfile(settings):
    return os.path.join(settings['outputfolder'], settings['name'], '%s.log' %settings['name'] )

loggeroperator = LoggerOperator()
    
def return_logger(name, include_thread):
    return loggeroperator.return_logger(name, include_thread)
    
def setup_loggers(settings):
    loggeroperator.startup(settings)

