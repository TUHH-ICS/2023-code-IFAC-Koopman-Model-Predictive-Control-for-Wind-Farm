function dataSOFWA = readdmdinformation(dirName)

dataSOFWA = cell(length(dirName),1);
for j = 1:length(dirName)
    
    [~, pitch] = readPitchData(strcat(dirName{j},'/1000/bladePitch')); %pitch=pitch{:};
    
    [~,time,~,~,powerGenerator] = readTurbineOutputGlobal(dirName{j},'generatorPower',inf);
    [~,~,~,~,rotSpeed] = readTurbineOutputGlobal(dirName{j},'rotorSpeed',inf);
    [~,~,~,~,rotorAzimuth] = readTurbineOutputGlobal(dirName{j},'rotorAzimuth',inf);
    [~,~,~,~,nacelleYaw] = readTurbineOutputGlobal(dirName{j}, 'nacelleYaw',inf);
    
    dataSOFWA{j}.time = time;
    dataSOFWA{j}.pitch = pitch;
    dataSOFWA{j}.powerGenerator = powerGenerator;
    dataSOFWA{j}.rotSpeed = rotSpeed;
    dataSOFWA{j}.rotorAzimuth = rotorAzimuth;
    dataSOFWA{j}.nacelleYaw = nacelleYaw;
end
