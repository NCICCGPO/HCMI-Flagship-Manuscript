import pandas as pd
from io import StringIO
import sys
pd.set_option('mode.use_inf_as_na', True)

def parse(log):
    table = []
    for line in log:
        if 'ProgressMeter' in line \
                and not 'raversal' in line:
            table.append(line)
    table = '\n'.join(table)
    return table


def clean(data):
    data = data[['Current', 'Locus',
       'Elapsed', 'Minutes', 'Regions', 'Processed', 'Regions/Minute']]
    data.columns = ['Current_Locus',
                    'Elapsed_Minutes', 'Regions_Processed',
                    'Regions/Minute', '1', '2', '3']
    return data[['Current_Locus',
       'Elapsed_Minutes', 'Regions_Processed', 'Regions/Minute']]


def add_details(progress):
    progress['prior_Regions_Processed'] = progress.Regions_Processed.shift(1)
    progress['prior_Elapsed_Minutes'] = progress.Elapsed_Minutes.shift(1)
    progress.loc[0, ['prior_Regions_Processed', 'prior_Elapsed_Minutes']] = 0
    progress['current_Regions_Processed'] = progress.Regions_Processed - progress['prior_Regions_Processed']
    progress['current_Elapsed_Minutes'] = progress.Elapsed_Minutes - progress['prior_Elapsed_Minutes']
    progress['regions_per_minute'] = progress['current_Regions_Processed'] / progress['current_Elapsed_Minutes']
    return progress

def merge(all):
    all['Regions_Processed'] = all['current_Regions_Processed'].cumsum()
    all['Elapsed_Minutes'] = all['current_Elapsed_Minutes'].cumsum()
    all['Regions/Minute'] = all['Regions_Processed'] / all['Elapsed_Minutes']
    return all

def get_progress(log_file):
    with open(log_file) as log:
        table = parse(log)
    data = pd.read_csv(StringIO(table), delimiter=r"\s+")
    return data

def parse_file(log_file):
    data = get_progress(log_file)
    data = clean(data)
    data = add_details(data)
    data.dropna(subset=['regions_per_minute'], inplace=True)
    return data


def main():
    log_files = sys.argv[1:-1]
    out_file = sys.argv[-1]
    dfs = []
    for log_file in log_files:
        data = parse_file(log_file)
        dfs.append(data)
    all = pd.concat(dfs, ignore_index=True)
    all = merge(all)
    all.to_csv(out_file, index=False)
    result = all.loc[all.shape[0]-1, 'Regions/Minute']
    sys.stdout.write(str(result))


if __name__ == "__main__":
    main()
