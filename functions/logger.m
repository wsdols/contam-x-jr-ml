classdef (Abstract) logger
    properties(Constant)
        LOG_RESULTS = 0x01;         % Log results
        LOG_PROJECT = 0x02;         % Log project definition
        LOG_SIM_AF = 0x04;          % Log airflow calculation
        LOG_SIM_AF_DETAILS = 0x08;  % Log details
        LOG_AFE = 0x10
    end
end