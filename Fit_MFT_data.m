function [x_fit, ci] = Fit_MFT_data(x0,Tarr,Barr,params,params_compare,params_indx,bounds_mft)
    
func_hand = @MFT_wrapper;

x0(5) = x0(5)/1e24;

options = optimset('MaxFunEvals',10000,'MaxIter',10000,'FunValCheck','on','Display','iter-detailed');
arr_struct.Tarr = Tarr;
arr_struct.Barr = Barr;
arr_struct.params_compare = params_compare;
arr_struct.params_indx = params_indx;
[x_fit,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(func_hand,x0,arr_struct,params,bounds_mft(1,:),bounds_mft(2,:),options);

ci = nlparci(x_fit,residual,'Jacobian',jacobian);

MFT_data = func_hand(x_fit,arr_struct);

figure
r = floor(sqrt(size(params,2)));
c = ceil(size(params,2)/r);
for k = 1:size(params,2)
    subplot(r,c,k)
    plot(Tarr,params(:,k),'kd');
    hold on
    plot(Tarr,MFT_data(:,k),'linestyle','-');
end                

disp(strcat(['Tc: ' num2str(x_fit(1)) ' ± ' num2str(x_fit(1)-ci(1,1))]));
disp(strcat(['gj: ' num2str(x_fit(2)) ' ± ' num2str(x_fit(2)-ci(2,1))]));
disp(strcat(['J: ' num2str(x_fit(3)) ' ± ' num2str(x_fit(3)-ci(3,1))]));
disp(strcat(['thetaD: ' num2str(x_fit(4)) ' ± ' num2str(x_fit(4)-ci(4,1))]));
disp(strcat(['Ns: ' num2str(x_fit(5)*1e24) ' ± ' num2str(x_fit(5)*1e24-ci(5,1)*1e24)]));
disp(strcat(['M: ' num2str(x_fit(6)) ' ± ' num2str(x_fit(6)-ci(6,1))]));
disp(strcat(['gamma_e: ' num2str(x_fit(7)) ' ± ' num2str(x_fit(7)-ci(7,1))]));

x_fit(5) = x_fit(5)*1e24;
end

function params = MFT_wrapper(x0,arr_struct)
    x0(5) = x0(5)*1e24;
    [Cp, DT, S, DS, mag]  = MFT_model(x0,arr_struct.Tarr,arr_struct.Barr);

    params = [];
    for i = 1:length(arr_struct.params_compare)
        switch arr_struct.params_compare{i}
            case 'Cp'
                data = Cp;
            case 'DT'
                data = DT;
            case 'S'
                data = S;
            case 'DS'
                data = DS;
            case 'mag'
                data = mag;
        end
        params = [params, data(arr_struct.params_indx(i),:)'];
    end
end