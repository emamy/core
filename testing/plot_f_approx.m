function plot_f_approx( model, model_data )
%PLOT_F_APPROX Summary of this function goes here
%   Detailed explanation goes here

    
    %TODO: continue here!
    xi = model_data.snapshots(:,:);
    fxi = model_data.f_values(:,:);
    
    vect = compile_triplevect(model, model_data);
    

end

