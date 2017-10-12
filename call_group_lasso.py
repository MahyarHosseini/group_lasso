import pandas
import numpy
import group_lasso


def read_file(file_addr, target_column_name, skip_column_list):
    df = pandas.read_excel(file_addr)
    columns = df.columns
    for col in df:
        if col not in skip_column_list:
            df[col] = [float(i) for i in df[col].values]
        else:
            if col == 'Exon_Malueka_Category':
                temp = []
                for val in df[col].values:
                    if val == 'A':
                        temp.append(1)
                    else:#val == 'C'
                        temp.append(0)
                df[col] = temp
          

    target = df[target_column_name].values
    del df[target_column_name]
    data = df.values#list of rows of the input file
 
    return [df.columns, data, target]


if __name__ == '__main__':
    from sklearn import datasets

    #Each columns of the data would be labeled with a number (columns with the same number are
    #in the same group)
    #Please label the columns after you printed the "df.columns" since the order might be
    #different from the input file
    #Skip to label the column that contains the target values
    #group_assigned = [0, 1, 1, 1, 1, 1, 2, 3, 4, 4, 5, 61, 7, 7, 7, 62, 63, 64, 65, 66, 67, 8, 9, 9, 9, 9, 9, 68, 69, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 3, 3, 3, 12, 12, 9, 13, 13, 13, 13, 13, 14]#Contains group label 
    group_assigned = [0, 1, 1, 1, 1, 1, 2, 3, 4, 4, 5, 6, 7, 7, 7, 62, 6, 6, 6, 6, 6, 8, 3, 3, 3, 3, 3, 6, 6, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 3, 3, 3, 12, 12, 9, 13, 13, 13, 13, 13, 14]#Contains group label 
#    group_assigned = [0, 1, 1, 1, 1, 1, 2, 3, 4, 4, 5, 6, 3, 3, 3, 6, 6, 6, 6, 6, 6, 3, 9, 9, 9, 9, 9, 6, 6, 10, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 12, 12, 9, 13, 13, 13, 13, 13, 14]#Contains group label 

    alpha = 0.15
    file_addr = '/home/seyedmah/Desktop/normalized_data_sep12.xlsx'
    target_column_name = 'skip_percentage'
    skip_column_list = ['Exon_Malueka_Category']
    columns, data, target = read_file(file_addr, target_column_name, skip_column_list)


    X = numpy.r_[data]
#    group_assigned = numpy.r_[numpy.arange(X.shape[1])]


    print '\nlen groups: ', len(group_assigned)
    print 'len columns (excluding target): ', len(columns)
    print 'Alpha = ', alpha
    print '\ngroup_assigned: ', group_assigned
   
    #X = numpy.r_[data]

    column_group_map = {}
    for i in range(len(group_assigned)):
        g = group_assigned[i]
        c = columns.values[i]
        if g in column_group_map:
            column_group_map[g].append(c)
        else:
            column_group_map[g] = [c]
    
    print '\nGroup Map:\n'
    for key in column_group_map:
        print key , ': ', column_group_map[key]
    print

    
    y = target
    groups = numpy.r_[group_assigned]
    coefs = group_lasso.group_lasso(X, y, alpha, groups, verbose=True, max_iter = 1000)
    print '\nKKT conditions verified:', group_lasso.check_kkt(X, y, coefs, alpha, groups), '\n'
   
    #Printing Non-zero Coefficents
    zero_flag = True
    for i in range(len(coefs)):
        if coefs[i] != 0:
            zero_flag = False
            print coefs[i] , ' * ', columns.values[i], '\t(Group =', group_assigned[i], ')'
    if zero_flag:
        print '\n All coefficents are zero!'

