def convert_csv_to_colored_xlsx(path_to_csv, data_col='b', from_col='a', to_col='e'):
    """
    :param path_to_csv: csv file to convert to a formatted xlsx
    :param number_of_loci: how many loci were analyzed
    :param data_col: column based on which the formatting will be done
    :param from_col: first column that should be formatted
    :param to_col: last column that should be formatted
    :return: creates a formatted xlsx file according to the
    """
    import pandas as pd

    df = pd.read_csv(path_to_csv)
    number_of_loci = df.shape[0]

    path_to_xlsx = path_to_csv.replace('csv', 'xlsx')
    writer = pd.ExcelWriter(path_to_xlsx, engine='xlsxwriter')
    df.to_excel(writer, sheet_name='Sheet1', index=None)

    workbook = writer.book
    worksheet = writer.sheets['Sheet1']

    # parameters adjustment
    number_of_loci += 1  # shift by 1 due to the header (start formatting from row number 2)
    data_col = data_col.upper()
    from_col = from_col.upper()
    to_col = to_col.upper()

    # table header formatting
    worksheet.conditional_format(f'{from_col}1:{to_col}1',
                                 {'type': 'formula',
                                  'criteria': '=$A1="locus"',
                                  'format': workbook.add_format({'bg_color': '#FFEB9C',
                                                                 'bold': True,
                                                                 'italic': True})})

    # category : minimal value, background color, font color
    categories = [[0.7, '#28A228', '#000000'],
                  [0.5, '#32C732', '#006100'],
                  [0.3, '#77DD77', '#006100'],
                  [0.1, '#B1D8B7', '#751100'],
                  [-1, '#FF9A8A', '#751100']]

    for category in categories:
        minimal_value, bg_color, font_color = category
        worksheet.conditional_format(f'{from_col}2:{to_col}{number_of_loci}',
                                     {'type': 'formula',
                                      'criteria': f'=${data_col}2>{minimal_value}',
                                      'format': workbook.add_format({'bg_color': bg_color,
                                                                     'font_color': font_color})})

    writer.save()

'''
if __name__ == '__main__':
    # test
    path_to_csv = 'concensus_predictions_full.csv'
    number_of_loci = 4305
    convert_csv_to_colored_xlsx(path_to_csv, number_of_loci)

'''