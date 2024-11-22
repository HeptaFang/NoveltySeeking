
ncutoff=15; % minimum number of cells to analyse session
bin=1; % bin of rasters used in GLM (in ms)
tau1=10; % synaptic integration time constant (in ms)
T1=4*tau1; % cutoff on the sum for the neuron total input
tau_all1=[0:T1];
tt_start=1+T1;
foldername='io_files_ovbins';

% Settings used previously to infer the Ising model parameters (used here as initial conditions to infer the GLM parameters):
bin_Ising=5; % bin of rasters (in ms)
regularization_Ising='L2'; 

%~~~~~~~~
ActiveNeurons_all=importdata(['io_files_ovbins/RasterNeurons_bin',num2str(bin),'.txt']);
ActiveBins_all=importdata(['io_files_ovbins/RasterTimeBins_bin',num2str(bin),'.txt']);
NB_all=importdata(['io_files_ovbins/N_B_bin',num2str(bin),'.txt']);
sessions_B=importdata(['io_files_ovbins/sessions_B_bin',num2str(bin),'.txt']);
load('RMcellist.mat')
%~~~~~~~~

% Extract nb. sessions and ind_sessions:
nrows	= size(infocell,2);
sessionNames = cell(nrows,1);
for i = 1:nrows
	sessionNames{i,1} = infocell(i).cellname(1:8);
end
sessions = unique(sessionNames,'stable');
nsessions = size(sessions,1);
ind_sessions=cell(nsessions,1);
for i=1:nsessions
    ind_sessions{i,1}=find(ismember(sessionNames,sessions{i,1}));
end

% Extract nb. cells and duration (in ms) for each session:
[ncells,raster_time_end]=nb_cells_times();
selected_sessions=find(ncells>=ncutoff);

options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'HessianFcn','objective','Display','iter');  

for sel=1:length(selected_sessions)
       
    session_name=sessions_B.textdata{sel,1}
    N=NB_all(sel,1);
    B=NB_all(sel,2);
    nb_par=N^2;
    
    ind_s=selected_sessions(sel,1);
    hJ_dhdJ=importdata(['io_files_ovbins/',sessions{ind_s,1},'_',regularization_Ising,'_bin',num2str(bin_Ising),'_learn.j']);
    J_Ising=zeros(N,N);
    cc=0;
    for i1=1:N-1
        for i2=i1+1:N
            cc=cc+1;
            J_Ising(i1,i2)=hJ_dhdJ(N+cc,1);
            J_Ising(i2,i1)=J_Ising(i1,i2);
        end
    end

    ActiveNeurons=ActiveNeurons_all(sel,:);
    ActiveNeurons=ActiveNeurons(~isnan(ActiveNeurons));
    ActiveBins=ActiveBins_all(sel,:);
    ActiveBins=ActiveBins(~isnan(ActiveBins));

    raster=zeros(B,N);
    for ii=1:length(ActiveNeurons)
        neur_ind=ActiveNeurons(1,ii);
        bin_ind=ActiveBins(1,ii);
        raster(bin_ind,neur_ind)=1;
    end

    % Initial condition:
    par_vector0=zeros(nb_par,1);  % par_vector=[{hi for all i}, {Jij for all j, then for all i}]
    par_vector0([1:N],1)=hJ_dhdJ([1:N],1);
    for i1=1:N
        Jnondiag=J_Ising(i1,:);
        Jnondiag(i1)=[];
        par_vector0(N+(i1-1)*(N-1)+[1:N-1],1)=Jnondiag;
    end
    
    predj=zeros(N,B);
    for tt=tt_start:B        
        for jj=1:N
            predj(jj,tt)=exp(-tau_all1./tau1)*raster(tt-tau_all1,jj); 
        end  
    end

    fun=@(par_vector)minuslogL_grad_hess_fun(par_vector,B,N,tt_start,raster,predj);

    [par_opt,cost_opt,exitflag,output,grad_minuslogL,hess_minuslogL] = fminunc(fun,par_vector0,options);
  
    stdErr_par=sqrt(diag(inv(hess_minuslogL)));
    
    save(['GLM_MLE_results/',session_name,'_results'],'par_opt','cost_opt','stdErr_par','grad_minuslogL','hess_minuslogL')
       
end


