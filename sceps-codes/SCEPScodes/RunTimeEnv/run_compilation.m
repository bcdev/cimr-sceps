

%=== compiling to produce the runtime executables

mcc -d /obs/cjimenez/Work/Tmp/RunTimeEnvRes -m -R  singleCompThread /obs/cjimenez/Work/Studies/SCEPS2OpenSF/SCEPSCodes/SceGenMod/Modules/GeoInputs_Extract.m
eval( ['!cp /obs/cjimenez/Work/Tmp/RunTimeEnvRes/GeoInputs_Extract .'])

mcc -d /obs/cjimenez/Work/Tmp/RunTimeEnvRes -m -R  singleCompThread /obs/cjimenez/Work/Studies/SCEPS2OpenSF/SCEPSCodes/SceGenMod/Modules/Forward_Model.m
eval( ['!cp /obs/cjimenez/Work/Tmp/RunTimeEnvRes/Forward_Model .'])




