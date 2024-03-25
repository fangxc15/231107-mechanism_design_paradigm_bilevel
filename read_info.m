function [Para,Num] = read_info(xlsname,sheetname, Setting)

    [data_raw,name_raw] = xlsread(['Data/',xlsname],sheetname);
    Num_column = find(strcmp(name_raw, 'Num'));
    type_column = find(strcmp(name_raw, 'type'));
    qmin_column = find(strcmp(name_raw, 'qmin'));
    qmax_column = find(strcmp(name_raw, 'qmax'));
    xmin_column = find(strcmp(name_raw, 'xmin'));
    xmax_column = find(strcmp(name_raw, 'xmax'));
    distri_column = find(strcmp(name_raw, 'distribution'));
    node_column = find(strcmp(name_raw, 'node'));
    company_column = find(strcmp(name_raw, 'company'));
    xstd_column = find(strcmp(name_raw, 'xstd'));
    xpara_column = find(strcmp(name_raw, 'xpara'));
    ifstrategic_column = find(strcmp(name_raw, 'ifstrategic'));
    
    Para.xmax = data_raw(:,xmax_column)';
    Para.xmin = data_raw(:,xmin_column)';
    Para.type = data_raw(:,type_column)';
    Para.qmax = data_raw(:,qmax_column)';
    Para.qmin = data_raw(:,qmin_column)';
    Para.distribution = data_raw(:,distri_column)';
    Para.node = data_raw(:,node_column)';
    Para.company = data_raw(:,company_column)';
    Para.xstd = data_raw(:,xstd_column)';
    Para.xpara = data_raw(:,xpara_column)';
%     Para.ifstrategic = data_raw(:,ifstrategic_column)';
    
    Para.xworse(Para.type == 1) = Para.xmin(Para.type == 1);
    Para.xworse(Para.type == 2) = Para.xmax(Para.type == 2);

    
    % 还有一个概率密度函数
    Num.I = max(data_raw(:,Num_column)');
    Num.C = max(Para.company);
    Num.P = Setting.NumP;
    [Para,Num] = process_Para(Para,Num,Setting);

end 