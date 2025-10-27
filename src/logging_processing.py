# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 15:29:02 2025

@author: u300737
"""
import logging

def setup_logging(log_filename):
    # Create a logger instance
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
    # Clear existing handlers to avoid duplicate logs
    if logger.hasHandlers():
        logger.handlers.clear()
    
    # Create a new file handler in write mode
    file_handler = logging.FileHandler(log_filename, mode='w')
    
    # Define a log message format with date, time, level and message
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    
    # Add the handler to the logger
    logger.addHandler(file_handler)
    
    return logger