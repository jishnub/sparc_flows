
import subprocess ,time, os, shutil,read_params,glob,fnmatch,sys
import numpy as np,pyfits
import matplotlib.pyplot as plt
from termcolor import colored
import re, datetime

################################################################################
# Check if compiled since last modification
################################################################################

def modification_date(filename):
    t = os.path.getmtime(filename)
    return datetime.datetime.fromtimestamp(t)

def check_if_built():
    with open("makefile","r") as f: l=f.readlines()
    exec_name = [x for x in l if x.lower().startswith("command1")][0].split("=")[-1].strip()
    sparc_build_time = modification_date(exec_name)
    for corefile in ["params.i",'driver.f90','integrals.f90','physics2d.f90','derivatives.f90',
                 'physics.f90','initialize.f90','all_modules.f90','pml.f90','step.f90',
                 'bspline90_22.f90','damping.f90','displacement.f90','traveltimes.f90',
                 'kernels.f90']:
        file_time = modification_date(corefile)
        if file_time > sparc_build_time:
            print("{} has been modified, please run make before running master.py".format(corefile))
            sys.exit(1)

check_if_built()

################################################################################

Rsun=695.9895 # Mm
z,c,rho = np.loadtxt(read_params.get_solarmodel(),usecols=[0,1,2],unpack=True); z=(z-1)*Rsun
c/=100 # m/s

#~ Get shape
nx = read_params.get_nx()
ny = 1
nz = read_params.get_nz()

Lx = read_params.get_xlength()
dx = Lx/nx

x = np.linspace(-Lx/2,Lx/2,nx,endpoint=False)

################################################################################

np.set_printoptions(precision=5)

datadir=read_params.get_directory()

iter_var="psi"

num_linesearches = 6

def grad_command(algo="bfgs",eps=[0.01*i for i in range(1,num_linesearches+1)]):
    eps_str = ' '.join([str(i) for i in eps])
    return "python grad.py algo="+algo+" "+str(eps_str)

id_text = read_params.parse_cmd_line_params(key="id")
opt_algo = read_params.parse_cmd_line_params(key="opt_algo",default="bfgs")

data_command = "qsub -N data_{} data_forward.sh".format(id_text)
full_command = "qsub -N full_{} full.sh".format(id_text)
ls_command = "qsub -N ls_{} linesearch.sh".format(id_text)


def copy_updated_model(num,var="psi"):
    updated_model=os.path.join(datadir,'update','test_'+var+'_'+str(num)+'.fits')
    model_for_next_iter=os.path.join(datadir,'model_'+var+'_ls00.fits')
    shutil.copyfile(updated_model,model_for_next_iter)

    try:
        updated_model=os.path.join(datadir,'update','test_'+var+'_'+str(num)+'_coeffs.npz')
        model_for_next_iter=os.path.join(datadir,'model_'+var+'_ls00_coeffs.npz')
        shutil.copyfile(updated_model,model_for_next_iter)
    except IOError: pass

if os.path.exists('compute_data'):
    os.remove('compute_data')

if os.path.exists('linesearch'):
    os.remove('linesearch')

if os.path.exists('running_full'):
    os.remove('running_full')

running_ls = False
running_full = False

# eps_fig = plt.figure()
model_fig = plt.figure()

def get_eps_around_minimum(prev1_eps,prev1_misfits):

    #~ Previous iteration, save the plot
    # eps_fig.clf()
    # eps_fig.gca().plot(prev1_eps,prev1_misfits,'bo',zorder=1)
    # plt.xlabel("Step size",fontsize=20)
    # plt.ylabel(r"misfit",fontsize=20)
    # plt.savefig('step_size_selection.eps')

    p=np.polyfit(prev1_eps,prev1_misfits,2)

    step = 1e-2
    if p[0]<0 and not (np.diff(prev1_misfits)<0).all():
        print(colored("inverted parabolic fit, can't estimate step size","red"))
        exit()
    elif p[0]<0 and (np.diff(prev1_misfits)<0).all():
        eps = np.array(prev1_eps)*10
    elif p[0]>0:
        step = - p[1]/(2*p[0])
        eps = np.array([step*(0.8+0.1*i) for i in range(num_linesearches)])
    elif p[0]==0:
        step = -p[2]/p[1]
        eps = np.array([step*(0.8+0.1*i) for i in range(num_linesearches)])


    pred_misfit = np.polyval(p,eps)
    print("Estimated step size",eps)
    print("Predicted misfit",pred_misfit)
    # plt.plot(eps,pred_misfit,'ro',zorder=1)

    #~ Plot the fit polynomial line, this should pass through the points and have a minimum
    #~ around the predicted points
    # epsfine = np.linspace(min(eps[0],prev1_eps[0]),max(eps[-1],prev1_eps[-1]),num=2000)
    # misfitfine = np.polyval(p,epsfine)
    # plt.plot(epsfine,misfitfine,'g-',zorder=0)

    # plt.savefig('step_size_selection.eps')
    # eps_fig.clf()

    return eps

def get_iter_status():
    ls_files=sorted(fnmatch.filter(os.listdir(os.path.join(datadir,"update")),'linesearch_[0-9][0-9]'))
    misfit_files=sorted(fnmatch.filter(os.listdir(os.path.join(datadir,"update")),'misfit_[0-9][0-9]'))
    ls_all_files=sorted(fnmatch.filter(os.listdir(os.path.join(datadir,"update")),'linesearch_all_[0-9][0-9]'))
    misfit_all_files=sorted(fnmatch.filter(os.listdir(os.path.join(datadir,"update")),'misfit_all_[0-9][0-9]'))
    iterno = len(misfit_files)-1
    return misfit_files,misfit_all_files,ls_files,ls_all_files,iterno

epslist_path=os.path.join(datadir,'epslist.npz')

def plot_true_and_model_psi():
    print("Plotting psi")
    model_fig.clf()
    model_fig.add_subplot(121)
    true_psi = np.squeeze(pyfits.getdata(read_params.get_true_psi_filename()))
    plt.pcolormesh(x,z,true_psi,cmap="RdBu_r")
    plt.xlim(-70,70)
    plt.ylim(-6,z[-1])
    plt.xlabel("x (Mm)",fontsize=16)
    plt.ylabel("z (Mm)",fontsize=16)
    plt.title("True",fontsize=14)
    plt.colorbar()
    model_fig.add_subplot(122)
    with np.load(os.path.join(datadir,"model_psi_ls00_coeffs.npz")) as f:
        back=f["back"].item()
    model_psi = np.squeeze(pyfits.getdata(os.path.join(datadir,"model_psi_ls00.fits")))
    plt.pcolormesh(x,z,model_psi-back,cmap="RdBu_r")
    plt.xlim(-70,70)
    plt.ylim(-6,z[-1])
    plt.xlabel("x (Mm)",fontsize=16)
    plt.ylabel("z (Mm)",fontsize=16)
    plt.title("Iter",fontsize=14)
    plt.colorbar()
    model_fig.set_size_inches(8,4)
    plt.tight_layout()
    plt.savefig(os.path.join(datadir,"model_psi_ls00.png"))

def print_status(all_lines):
    timesteps_regex = re.compile(r"NUMBER OF TIMESTEPS =\s+\d+")
    timesteps_line = timesteps_regex.findall(all_lines)
    if not timesteps_line: return (0,0,0,0)
    number_of_timesteps = int(timesteps_line[0].split()[-1])
    status_regex = re.compile(r"Iteration #,  Cpu Time and Real Time:\n\s+\d+\s \d+\.\d+")
    iter_lines = status_regex.findall(all_lines)
    if not iter_lines: return (0,number_of_timesteps,0,0)
    last_iter_line = iter_lines[-1]
    steps_done = int(last_iter_line.split()[-2])
    time_per_step = float(last_iter_line.split()[-1])
    steps_remaining = number_of_timesteps - steps_done
    estimated_time = round(steps_remaining*time_per_step/60)
    return (steps_done,number_of_timesteps,steps_done/number_of_timesteps*100,
    estimated_time)

for query in range(100000):

    #~ print("Query",query)
    qstat=subprocess.check_output(["qstat"]).split()

    job_id_indices = [qstat.index(i) for i in [x for x in qstat if x.endswith("daahpc1") or x.endswith("daahpc2")]]
    job_names = [qstat[i+1] for i in job_id_indices]
    job_statuses = [qstat[i+4] for i in job_id_indices]

    job_running = False
    for job_name,job_status in zip(job_names,job_statuses):
        if "_".join(job_name.split("_")[1:])==id_text and (job_status=="R" or job_status=="Q"):
            job_running = True
            what_is_running = job_name.split("_")[0]
            break

    if job_running:
        fmtstr = "{:10s}: {:5d} of {:<5d},{:3.0f}%, ETR {:3.0f} mins"
        if what_is_running=="data":
            outfile = os.path.join(datadir,"forward_src01_ls00","out_data_forward")
            if os.path.exists(outfile):
                with open(outfile,"r") as f: all_lines = f.read()
                print(fmtstr.format("Data",*print_status(all_lines)),end="\r")

        elif what_is_running=="full":
            outfile = os.path.join(datadir,"forward_src01_ls00","out_forward")
            if os.path.exists(outfile):
                with open(outfile,"r") as f: all_lines = f.read()
                steps_done,number_of_timesteps = print_status(all_lines)[:2]
                if steps_done<number_of_timesteps and number_of_timesteps>0:
                    print(fmtstr.format("Forward",*print_status(all_lines)),end="\r")
                elif steps_done==number_of_timesteps and number_of_timesteps>0:
                    outfile = os.path.join(datadir,"adjoint_src01","out_adjoint")
                    if os.path.exists(outfile):
                        with open(outfile,"r") as f: all_lines = f.read()
                        print(fmtstr.format("Adjoint",*print_status(all_lines)),end="\r")

        elif what_is_running=="ls":
            outfile = os.path.join(datadir,"forward_src01_ls01","out_linesearch")
            if os.path.exists(outfile):
                with open(outfile,"r") as f: all_lines = f.read()
                print(fmtstr.format("Linesearch",*print_status(all_lines)),end="\r")

        sys.stdout.flush()
        time.sleep(20)
        continue
    else: print()

    if not os.path.exists(os.path.join(datadir,"data","01.fits")):
        #~ no iterations done
        print(colored("Running data_forward","blue"))
        status=subprocess.call(data_command.split())
        assert status==0,"Error in running data forward"
        time.sleep(60)
        continue

    misfit_files,misfit_all_files,ls_files,ls_all_files,iterno = get_iter_status()
    #~ print(misfit_files,ls_files)
    print(colored("Iteration "+str(iterno),"green"))

    if iterno==-1:
        #~ Start of iterations
        plot_true_and_model_psi()
        print("#"*80)
        print(colored("Starting with iterations, running full for the first time","blue"))
        status=subprocess.call(full_command.split())
        assert status==0,"Error in running full"
        time.sleep(20)
        continue

    elif len(misfit_files)>len(ls_files):

        #~ Kernel computation should be over for this iteration
        #~ We should proceed to computing update and running linesearch

        plot_true_and_model_psi()

        running_full = False
        if running_ls:
            print(colored("It seems the linesearch file wasn't generated. Check if previous linesearch finished correctly","red"))
            exit()
        misfit_from_prev_iter = np.sum(np.loadtxt(os.path.join(datadir,"update",
                                misfit_files[-1]),usecols=[2]))
        print("Misfit from previous iteration",misfit_from_prev_iter)

        # Check if misfit has changed substantially
        misfit_from_two_iters_back = 0
        if len(misfit_files)>1:
            misfit_from_two_iters_back = np.sum(np.loadtxt(os.path.join(
                                            datadir,"update",misfit_files[-2]),
                                            usecols=[2]))
        change_in_misfit = misfit_from_prev_iter - misfit_from_two_iters_back
        if len(misfit_files)>1 and change_in_misfit>0:
            print("travel-time misfit increasing, stopping")
            break
        if abs(change_in_misfit)/misfit_from_prev_iter<1e-3:
            print("Misfit from 2 iters back {}".format(misfit_from_two_iters_back))
            print("Misfit this iter {}".format(misfit_from_prev_iter))
            print("Relative change in misfit {}".format(abs(change_in_misfit)/misfit_from_prev_iter))
            print("Change in travel-time misfit in very small, stopping")
            break


        #~ Need to run linesearch for this iteration
        #~ check if eps for iteration exists
        try:
            epslist=np.load(epslist_path)
            iters_list = epslist.files

            if str(iterno) in iters_list:
                #~ If it exists, this means linesearch files were removed
                #~ This probably means that the linesearch misfits were monotonic

                ls_details = epslist[str(iterno)]

                #~ print(ls_details)

                prev1_done = ls_details[1].all()
                prev2_done = ls_details[3].all()

                if prev1_done:
                    print("Computing step size from previous linesearch")
                    eps_prev = ls_details[0]
                    lsmisfit_prev =  ls_details[1]
                    eps = get_eps_around_minimum(eps_prev,lsmisfit_prev)

                elif not(prev1_done) and prev2_done:
                    print("Computing step size form 2nd last linesearch")
                    eps_prev = ls_details[2]
                    lsmisfit_prev =  ls_details[3]
                    eps = get_eps_around_minimum(eps_prev,lsmisfit_prev)

                elif (not prev1_done) and (not prev2_done):
                    #~ Probably interrupted, so misfits not saved but step sizes might be recorded
                    #~ Restart with one set of step sizes
                    if ls_details[0].any():
                        print("Reusing previous step size")
                        eps = ls_details[0]
                    elif ls_details[2].any():
                        print("Reusing step size from 2nd last iteration")
                        eps = ls_details[2]
                    else:
                        print("No step size information recorded",)
                        prev_iternos= np.array([int(i) for i in iters_list if i!=str(iterno)])
                        if len(prev_iternos)!=0:
                            closest_iter = str(prev_iternos[abs(prev_iternos-iterno).argmin()])
                            eps=epslist[closest_iter][0]
                            print("using step size from iteration",closest_iter)
                        else:
                            print("using arbitary step sizes")
                            eps=[1e-2*i for i in range(1,num_linesearches+1)]
            else:
                print("Linesearch for this iteration not done yet, choosing step size from iteration",)
                #~ Iteration not yet done
                #~ Choose the value corresponding to the closest iteration
                prev_iternos= np.array([int(i) for i in iters_list])
                closest_iter = str(prev_iternos[abs(prev_iternos-iterno).argmin()])
                print(closest_iter)
                eps=epslist[closest_iter][0]

        except IOError:
            #~ If the file doesn't exist, no linesearches have been carried out
            #~ Assign arbitrary small value
            print("Using arbitary step sizes")
            eps=np.array([1e-2*i for i in range(1,num_linesearches+1)])
        #~ Run grad
        gradcmd = grad_command(algo=opt_algo,eps=eps)
        def exp_notation(a):
            b=a.split()
            for index,word in enumerate(b):
                try:
                    b[index] = "{:2.1E}".format(float(word))
                except ValueError: pass
            return ' '.join(b)

        print("Running",exp_notation(gradcmd))
        status=subprocess.call(gradcmd.split())
        assert status==0,"Error in running grad"

        #~ Run linesearch
        print(colored("Running linesearch","blue"))
        status=subprocess.call(ls_command.split())
        assert status==0,"Error in running linesearch"
        running_ls = True
        time.sleep(60)
        continue

    elif iterno>=0 and len(misfit_files)==len(ls_files):

        #~ Linesearch should be over, and we can do two things -
        #~ (a) run linesearch again with different step size if misfit minimum is not found
        #~ (b) run full if misfit minimum is found

        running_ls = False
        if running_full:
            print(colored("It seems the full misfit was not generated. Check if the previous full.sh finished without errors","red"))
            exit()
        #~ check linesearch misfits
        ls_latest = os.path.join(datadir,"update",ls_files[-1])
        # ls_all_latest = os.path.join(datadir,"update",ls_all_files[-1])

        try:
            lsdata = np.loadtxt(ls_latest)
        except IOError:
            print(colored("Could not load "+ls_latest+", check if file exists","red"))
            exit()

        assert len(np.where(lsdata[:,0]==1)[0])==num_linesearches,"Not all linesearches finished correctly"
        nmasterpixels = lsdata.shape[0]//num_linesearches
        misfit=np.array([sum(lsdata[i*nmasterpixels:(i+1)*nmasterpixels,2]) for i in range(num_linesearches)])

        print("Linesearch misfits",misfit)
        if np.isnan(misfit).any():
            print("Nan found, quitting")
            exit()
        #~ Check for misfit minimum
        misfit_diff = np.diff(misfit)
        misfit_diff_sign_change_locs = np.where(np.diff(np.sign(misfit_diff)))[0]
        misfit_diff_changes_sign = len(misfit_diff_sign_change_locs)>0
        if misfit_diff_changes_sign:
            misfit_diff_sign_change_ind = misfit_diff_sign_change_locs[0]
            # Ensure that it is a simple minimum
            if not ((np.sign(misfit_diff)[:misfit_diff_sign_change_ind+1]==-1).all() and
                    (np.sign(misfit_diff)[misfit_diff_sign_change_ind+1:]==1).all()):
                    print("Not a simple minimum, quitting")
                    exit()

        #~ Update epslist
        epslist=np.load(epslist_path)
        if str(iterno) in epslist.files:
            ls_details = epslist[str(iterno)]
            ls_details[1] = misfit
        else:
            ls_details = np.zeros_like(epslist[epslist.files[0]])
        epslist = {i:epslist[i] for i in epslist.files}
        epslist.update({str(iterno):ls_details})
        #~ print(epslist[str(iterno)])
        np.savez(epslist_path,**epslist)

        #~ Copy best linesearch, or remove prev linesearch file
        if misfit_diff_changes_sign:
            #~ Get index of minimum. Add two ones to the minimum of the diff array
            #~ One because the sign change occurs at the index before the minimum
            #~ The other because the test model counting starts from 1 instead of 0
            min_ls_index = misfit_diff_sign_change_locs[0]+1+1
            print("Misfit sign change detected")
            print("model #"+str(min_ls_index)+" has minimum misfit")
            copy_updated_model(min_ls_index,var=iter_var)

            plt.clf() # Refresh the stepsize-misfit plot

            #~ Run full.sh
            print("#"*80)
            plot_true_and_model_psi()
            print(colored("Running full","blue"))
            status=subprocess.call(full_command.split())
            running_full=True
            assert status==0,"Error in running full"
            time.sleep(60)
        else:
            if misfit_diff[0]>0: print("Linesearch strictly increasing")
            elif misfit_diff[0]<0: print("Linesearch strictly decreasing")

            #~ Monotonic linesearches don't count, remove previous linesearch file
            print("Removing linesearch file")

            ls_latest_renamed=ls_latest+".rnm"
            os.rename(ls_latest,ls_latest_renamed)

            # ls_all_latest_renamed= ls_all_latest+".rnm"
            # os.rename(ls_all_latest,ls_all_latest_renamed)

        continue
