function spot_struct = truncateParticleFields(spot_struct,has3DSpotInfo,hasAPInfo)

    spot_struct = spot_struct(~isnan([spot_struct.particleID]));
    spot_struct = rmfield(spot_struct,'spotFrames');
    
    for i = 1:length(spot_struct)
      fluo = spot_struct(i).fluo;
      nanFilter = ~isnan(fluo);
      
      spot_struct(i).fluo = fluo(nanFilter);
      spot_struct(i).time = spot_struct(i).time(nanFilter);
      spot_struct(i).frames = spot_struct(i).frames(nanFilter);
      
      spot_struct(i).xPosParticle = spot_struct(i).xPosParticle(nanFilter);
      spot_struct(i).yPosParticle = spot_struct(i).yPosParticle(nanFilter);
      spot_struct(i).zPosParticle = spot_struct(i).zPosParticle(nanFilter);      
      
      spot_struct(i).xPosNucleus = spot_struct(i).xPosNucleus(nanFilter);
      spot_struct(i).yPosNucleus = spot_struct(i).yPosNucleus(nanFilter);
      spot_struct(i).rawNCPprotein = spot_struct(i).rawNCPprotein(nanFilter);
      
      spot_struct(i).FrameQCFlags = spot_struct(i).FrameQCFlags(nanFilter);
      
      if has3DSpotInfo
        spot_struct(i).xPosParticle3D = spot_struct(i).xPosParticle3D(nanFilter);
        spot_struct(i).yPosParticle3D = spot_struct(i).yPosParticle3D(nanFilter);
        spot_struct(i).zPosParticle3D = spot_struct(i).zPosParticle3D(nanFilter); 
%         spot_struct(i).zPosParticleDetrended3D = spot_struct(i).zPosParticleDetrended3D(nanFilter);
        spot_struct(i).fluo3D = spot_struct(i).fluo3D(nanFilter);
      end
      if hasAPInfo 
        spot_struct(i).APPosParticle = spot_struct(i).APPosParticle(nanFilter);
      end

      spot_struct(i).zPosParticleDetrended = spot_struct(i).zPosParticleDetrended(nanFilter);
      
      spot_struct(i).fluoOffset = spot_struct(i).fluoOffset(nanFilter);       
    end