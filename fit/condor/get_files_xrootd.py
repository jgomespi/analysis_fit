import subprocess, os, re


def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def generate_path(dataset,  crab_folder, n_folders):

    cmds = [ 
    f'xrdfs xrootd-redir.ultralight.org ls -u /store/group/uerj/mabarros/Charmonium/{dataset}_AOD/{crab_folder}/{i:04}/' for i in range(n_folders)
    ]
    #print(cmds)

    cat = ""
    out_file = dataset + "_path.txt"
    for i in cmds:
        os.system(f"{i} > {i.split('/')[-2]}")
        cat += f" {i.split('/')[-2]}"

    os.system(f"cat {cat} > {out_file}")
    os.system(f"rm -rf {cat}")
    file_list = open(out_file, 'r').readlines()

    for idx, f in enumerate(file_list):
        file_list[idx] = re.sub('transfer-\d*', 'xrootd-redir', f)

    list(set(file_list))
    file_list.sort(key=natural_keys)

    final_list = []
    for i in file_list:
        if i not in final_list:
            final_list.append(i)

    with open(out_file, 'w') as f:
        for i in final_list:
            f.write(i)

if __name__ == '__main__':


    dataset = ['CharmoniumRun2016B_21Feb2020_ver2_UL2016_HIPM-v1',
               'CharmoniumRun2016C_21Feb2020_UL2016_HIPM-v1', 
               'CharmoniumRun2016D_21Feb2020_UL2016_HIPM-v1',
               'CharmoniumRun2016E_21Feb2020_UL2016_HIPM-v1',
               'CharmoniumRun2016F_21Feb2020_UL2016_HIPM-v1',
               'CharmoniumRun2016F_21Feb2020_UL2016-v1',
               'CharmoniumRun2016G_21Feb2020_UL2016-v1',
               'CharmoniumRun2016H_21Feb2020_UL2016-v1']
    crab_folder = ['221009_034842', '221009_021601', '221009_034847', '221009_034852',
                   '221009_034858', '221009_034903', '221009_034908', '221009_040921']
 
    n_folders = [4, 2, 3, 3, 2, 1, 5, 6]
    
    for d, c, n in zip(dataset, crab_folder, n_folders):
        generate_path(d, c, n)     
