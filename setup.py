def create_directories(datadir,num_src,num_ls):

    for src in range(1,num_src+1):

        for job in range(num_ls+1):
        
            (datadir/"forward_src{:02d}_ls{:02d}".format(src,job)).mkdir(parents=True,exist_ok=True)

        (datadir/"adjoint_src{:02d}".format(src)).mkdir(parents=True,exist_ok=True)

    for subdir in ["kernel","update","status","tt"]:
        (datadir/subdir).mkdir(parents=True,exist_ok=True)
