function compileDKE
%   Compile DKE
%   Just run function to compile DKE into an executable for Mac

mkdir('DKE_compiled');
mcc -m -v  dke.m -a ./mac_mex -I ./mac_mex -d ./DKE_compiled -R '-logfile,dke_log.log'
end