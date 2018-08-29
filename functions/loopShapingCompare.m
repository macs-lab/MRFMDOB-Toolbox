function [allmargin_wComp,allmargin_woComp] = ...
    loopShapingCompare(P,C_baseline,C_wComp,Ts,flag_save_plot,Sfig_num,SW_Sonly)
% function [allmargin_wComp,allmargin_woComp] = ...
%     loopShapingCompare(P,C_baseline,C_wComp,Ts,flag_save_plot)
% * Evaluates frequency-domain loop shaping performances
% * Compares two servo designs: 
%   - one with a baseline controller C_baseline
%   - the second, C_wComp, with some additional compensation algorithms.
%
% ============================================================
%   Copyright (c) 2008-, Xu Chen, xchen@engr.uconn.edu 
%   Author(s): Xu Chen
%   initial version: 2015-01-04
% ============================================================
if nargin < 7
    SW_Sonly = 0;
end
if nargin < 6
    Sfig_num = [];
end
if nargin < 5
    flag_save_plot = 'nosave';
end

if length(fieldnames(C_baseline)) ~= 3 
    % a SISO design
    C_DS.Cvcm = C_baseline;
    C_DS.Cma = tf(0,1,Ts);
    C_DS.P2hat = tf(0,1,Ts);
else % a dual-input single-output (DISO) design
    C_DS = C_baseline; 
end
if strcmp(class(P),'frd')
    SW_FRD_DATA = 1;
else
    SW_FRD_DATA = 0;
end
%% Set up
SW_TUNE                     = 0;
SW_MANUAL_BODE              = 1;
z                           = tf('z',Ts);

if size(P,2) == 2
    FLAG_DSA = 1;
    C_woComp = [series(1+series(C_DS.Cma,C_DS.P2hat),C_DS.Cvcm);C_DS.Cma];
else
    FLAG_DSA = 0;
    C_woComp = C_DS.Cvcm;
end

bode_opt                = bodeoptions;
bode_opt.Xlim           = [50 0.5/Ts];
bode_opt.FreqUnits      = 'Hz';
bode_opt.PhaseWrapping  = 'On';
bode_opt.Grid           = 'On';
%% Transfer function computation
L_wComp          = P*C_wComp;
%                               loop tranfer func
S_wComp          = 1/(L_wComp + 1);
%                               sensitivity func
L_woComp         = P*C_woComp;
%                               loop tranfer func
S_woComp         = 1/(L_woComp+1);
%                               sensitivity func

% MATLAB's own bode function does not calculate L_wComp right
% need to do a manual bode

if SW_FRD_DATA
    w = P.Frequency;
else
    %
    % get a fine frequency grid using L_wComp
    %
%     [~,~,w]          = bode(L_wComp);
    [~,~,w]          = bode(C_wComp(1)); % 2015-01-10
    if length(w)<200
        w = linspace(min(w),max(w),10000);
    end
end
if min(w/2/pi)<10
    xcale = [10,max(w/2/pi)];
else
    xcale = [min(w/2/pi),max(w/2/pi)];
end
%
% freqresp: loop TFs w/o compensator
if SW_FRD_DATA
    mag_L = abs(L_woComp.ResponseData);
    %     ph_L = phase(L_woComp.ResponseData);
else
    [mag_L,ph_L]            = bode(L_woComp,w);
end
freqresp_L              = freqresp(L_woComp,w);
freqresp_L              = freqresp_L(:);
freqresp_L_woComp        = freqresp_L;
freqresp_S_woComp        = 1./(1+freqresp_L);
%
% freqresp: loop TFs w/ compensator
% frequency response from matlab's bode, this is wrong, we won't use it
% we will take only the frequency vector
if ~SW_FRD_DATA
    [mag_L_wComp,ph_L_wComp]  = bode(L_wComp,w);
    mag_L                   = mag_L(:);
    mag_L_wComp              = mag_L_wComp(:);
    
    ph_L                    = ph_L(:);
    ph_L_wComp               = ph_L_wComp(:);
    %
    % freqresp: P
    freqresp_P         = freqresp(P,w);
    freqresp_P1         = freqresp_P(1,1,:);
    freqresp_P1         = freqresp_P1(:);
    try % if it is a DISO plant
        freqresp_P2         = freqresp_P(1,2,:);
        freqresp_P2        = freqresp_P2(:);
    catch
    end
end
%
% freqresp: C
freqresp_C_wComp         = freqresp(C_wComp,w);
freqresp_C1_wComp        = freqresp_C_wComp(1,1,:);
freqresp_C1_wComp        = freqresp_C1_wComp(:);
try % if it is a DISO plant
    freqresp_C2_wComp        = freqresp_C_wComp(2,1,:);
    freqresp_C2_wComp        = freqresp_C2_wComp(:);
catch
end
%
% manual freqresp
if SW_FRD_DATA
    freqresp_LComp_manual = freqresp(L_wComp,w);
else
    if FLAG_DSA
        freqresp_LComp_manual = freqresp_P1.*freqresp_C1_wComp + ...
            freqresp_P2.*freqresp_C2_wComp;
    else
        freqresp_LComp_manual = freqresp_P1.*freqresp_C1_wComp;
    end
end
mag_LComp_manual         = abs(freqresp_LComp_manual);
if ~SW_FRD_DATA
    ph_LComp_manual          = phase(freqresp_LComp_manual)/pi*180;
end
freqresp_L_wComp         = freqresp_LComp_manual;
freqresp_S_wComp         = 1./(1+freqresp_L_wComp);
%
% for the phase, transform to deg and then normalize to [-180,180]
if ~SW_FRD_DATA
    ph_L_woComp_manual       = phase(freqresp_L_woComp)*180/pi;
    ph_L_woComp_manual       = mod(ph_L_woComp_manual+180,360)-180;
    ph_L_wComp_manual        = phase(freqresp_L_wComp)*180/pi;
    ph_L_wComp_manual        = mod(ph_L_wComp_manual+180,360)-180;
end
%% Plot
if SW_TUNE
    %
    % check that matlab's bode gives the wrong response
    figure, subplot(211)
    semilogx(w/2/pi,mag2db(mag_LComp_manual))
    hold on
    semilogx(w/2/pi,mag2db(mag_L_wComp),'r--')
    xlim([min(w/2/pi),max(w/2/pi)])
    legend('manual','MATLAB bode','location','best')
    ylabel ('Magnitude (dB)')
    subplot(212)
    semilogx(w/2/pi,ph_LComp_manual)
    hold on
    semilogx(w/2/pi,ph_L_wComp,'r--')
    xlim([min(w/2/pi),max(w/2/pi)])
    xlabel ('Frequency (Hz)')
    ylabel ('Phase (deg)')
end
%
% sensitivity compare
if ~SW_FRD_DATA
    if isempty(Sfig_num)
        figure
    else
        figure(Sfig_num)
    end
    if SW_MANUAL_BODE
        semilogx(w/2/pi,mag2db(abs(freqresp_S_woComp)),'linewidth',2)
        hold on
        semilogx(w/2/pi,mag2db(abs(freqresp_S_wComp)),'r--','linewidth',2)
        xlim(xcale)
        xlabel ('Frequency (Hz)')
        ylabel ('Magnitude (dB)')
    else
        h =         bodeplot(S_woComp,'r:',bode_opt);
        setoptions(h,'PhaseVisible','Off');
        h =         findobj(gcf,'type','line');
        set         (h,'linewidth',2);
        hold on;
        bodeplot    (S_wComp);
        h =         findobj(gcf,'type','line');
        set         (h,'linewidth',2);
        hold off;
        h =         findobj(gcf,'type','axes');
        set         (h,'fontsize',8);
        pause       (.5);
    end
    legend      (...
        'w/o compensator',...
        'w/ compensator',...
        'Location', 'southwest')
    title       'Frequency response of the sensitivity function'
    grid off
else 
    if isempty(Sfig_num)
        figure
    else
        figure(Sfig_num)
    end
    sw_temp = 1;
    if 1
        xbodemag({S_wComp,S_woComp},2)
    else
        xbodeplot({S_wComp,S_woComp},2)
    end
    legend      (...
        'w/ compensator',...
        'w/o compensator',...
        'Location', 'Best')
    %     title       'Frequency response of the sensitivity function'
    grid off
    if ~sw_temp
        subplot(212)
        grid off
        subplot(211)
    end
    
end
if strncmpi(flag_save_plot,'save',4)
    hgsave('bode_sensitivity_compare')
    saveas(gcf,'bode_sensitivity_compare.emf')
    try
        print -painters -depsc bode_sensitivity_compare
    catch
    end
end
%
% nyquist plot of loop tfs, sometimes difficult to tell the change due to
% bad scaling
if ~SW_Sonly
    figure
    nyquist     (L_woComp, L_wComp, 'r--')
    legend      ('w/o compensator','w/ compensator','location','best')
    title       'Loop transfer function'
    if strncmpi(flag_save_plot,'save',4)
        hgsave('nyquist_compare')
        saveas(gcf,'nyquist_compare.emf')
        try
            print -painters -depsc nyquist_compare
        catch
        end
    end
    
    % figure
    % nyquist     (L_woComp, L_wComp, 'r--')
    % zoom on
    % legend      ('w/o compensator','w/ compensator','location','best')
    % title       'Loop transfer function'
    %
    % figure
    % nyquist     (L_woComp, L_wComp, 'r--')
    % zoom on
    % legend      ('w/o compensator','w/ compensator','location','best')
    % title       'Loop transfer function'
    %
    % equivalent cascaded feedback controller
    % figure, bodeplot(G_FB_eqval,bode_opt)
    % set_bode_line(1.5)
    % title 'Equivalent cascaded feedback controller'
    % if strncmpi(flag_save_plot,'save',4)
    %     hgsave('bode_equivalent_fb_addon_controller')
    %     saveas(gcf,'bode_equivalent_fb_addon_controller.emf')
    %     try
    %         print -painters -depsc bode_equivalent_fb_addon_controller
    %     catch
    %     end
    % end
    %
    % loop transfer functions
    if ~SW_FRD_DATA
        if SW_MANUAL_BODE
            figure,     subplot(211)
            semilogx(w/2/pi,mag2db(abs(freqresp_L_woComp)),'linewidth',2)
            hold on
            semilogx(w/2/pi,mag2db(abs(freqresp_L_wComp)),'r--','linewidth',2)
            legend      ('w/o compensator','w/ compensator','location','best')
            xlim(xcale)
            ylabel      ('Magnitude (dB)')
            subplot(212)
            semilogx    (w/2/pi,ph_L_woComp_manual,'linewidth',2)
            hold on
            semilogx    (w/2/pi,ph_L_wComp_manual,'r--','linewidth',2)
            xlim        (xcale)
            xlabel      ('Frequency (Hz)')
            ylabel      ('Phase (deg)')
        else
            figure,     bodeplot(L_woComp, L_wComp, 'r--', bode_opt)
            h =         findobj(gcf,'type','line');
            set         (h,'linewidth',1.5);
            legend      ('w/ compensator','w/o compensator','location','best')
            grid off
        end
    else
        figure,     xbodeplot({L_woComp, L_wComp})
        legend      ('w/ compensator','w/o compensator','location','best')
        title       'Loop transfer function'
    end
    if strncmpi(flag_save_plot,'save',4)
        hgsave('bode_loop_TF_compare')
        saveas(gcf,'bode_loop_TF_compare.emf')
        try
            print -painters -depsc bode_loop_TF_compare
        catch
        end
    end
end
%% Stability margins
if ~SW_FRD_DATA
    [numcl,dencl]           = tfdata(S_woComp,'v');
    [clzero1,clpole1,clk]   = tf2zp(numcl,dencl);
    % disp '==================w/o Comp===================='
    % disp('magnitude of closed loop poles)')
    % abs(clpole1)'
    
    [numcl,dencl]           = tfdata(S_wComp,'v');
    [clzero1,clpole1,clk]   = tf2zp(numcl,dencl);
    % disp '==================w Comp======================'
    % disp('magnitude of closed loop poles)')
    % abs(clpole1)'
    if 0
        figure,     pzplot(S_woComp, S_wComp, 'r--')
        legend      (...
            'w/o Comp',...
            'w/ Comp',...
            'location','best')
    end
    
    %
    % % Closed loop robust stability test
    % deltaG = ...
    %     ( L - z^(-d_DOB)/Ln_inv )/...
    %     ( z^(-d_DOB)/Ln_inv );
    % if 1
    %     if SW_MANUAL_BODE
    %         [~,~,w_ana] = bode((1+L)/(z^-d_DOB*deltaG));
    %         [~,~,w_Q] = bode(Q_PF);
    %         % get a fine frequency grid first
    %         w_append = [w_ana(:);w_Q(:)];
    %         w_append = sort(w_append);
    %         w_append = w_append(w_append>=10*2*pi);
    %
    %         [mag_robuAna,~] = bode((1+L)/(z^-d_DOB*deltaG),w_append);
    %         [mag_Q,~] = bode(Q_PF,w_append);
    %
    %         mag_robuAna = mag_robuAna(:);
    %         mag_Q = mag_Q(:);
    %         figure
    %         semilogx(w_append/2/pi,mag2db(mag_robuAna),'linewidth',2)
    %         hold on
    %         semilogx(w_append/2/pi,mag2db(mag_Q),'r--','linewidth',2)
    %         xlim([min(w_append/2/pi),max(w_append/2/pi)])
    %         xlabel('Frequency (Hz)')
    %         ylabel('Magnitude (dB)')
    %         grid
    %     else
    %         figure
    %         h=bodeplot((1+L)/(z^-d_DOB*deltaG),bode_opt);
    %         setoptions(h,'FreqUnits','Hz','PhaseVisible','Off');
    %         hold on;
    %         bodeplot(Q_PF,'r--');
    %         xlim([10,0.5/Ts])
    %         h = findobj(gcf,'type','line');
    %         set(h,'linewidth',2);
    %         h = findobj(gcf,'type','axes');
    %         set(h,'fontsize',8);
    %         hold off;grid on;
    %         pause(0.5)
    %     end
    % else
    %     figure
    %     xbodemag({(1+L)/(z^-d_DOB*deltaG),Q_PF});
    % end
    % legend('(1+L)/\Delta','Q','Location','best');
    % if strncmpi(flag_save_plot,'save',4)
    %     hgsave('bode_robust_stability')
    %     saveas(gcf,'bode_robust_stability.emf')
    %     try
    %         print -painters -depsc bode_robust_stability
    %     catch
    %     end
    % end
    allmargin_wComp = allmargin(L_wComp);
    allmargin_woComp = allmargin(L_woComp);
    disp('===============Stability Margins w/ compensator===============')
    allmargin_wComp.GMFrequency = allmargin_wComp.GMFrequency./2./pi;
    allmargin_wComp.PMFrequency = allmargin_wComp.PMFrequency./2./pi;
    allmargin_wComp.DMFrequency = allmargin_wComp.DMFrequency./2./pi;
    allmargin_woComp.GMFrequency = allmargin_woComp.GMFrequency./2./pi;
    allmargin_woComp.PMFrequency = allmargin_woComp.PMFrequency./2./pi;
    allmargin_woComp.DMFrequency = allmargin_woComp.DMFrequency./2./pi;
    
    disp(allmargin_wComp)
    disp('GainMargin:')
    disp(allmargin_wComp.GainMargin)
    disp('GMFrequency:')
    disp(allmargin_wComp.GMFrequency)
    disp('PhaseMargin:')
    disp(allmargin_wComp.PhaseMargin)
    disp('PMFrequency:')
    disp(allmargin_wComp.PMFrequency)
    disp('===============Stability Margins w/o compensator===============')
    disp(allmargin_woComp)
    disp('GainMargin:')
    disp(allmargin_woComp.GainMargin)
    disp('GMFrequency:')
    disp(allmargin_woComp.GMFrequency)
    disp('PhaseMargin:')
    disp(allmargin_woComp.PhaseMargin)
    disp('PMFrequency:')
    disp(allmargin_woComp.PMFrequency)
    disp('===============================================')
else
    allmargin_wComp = [];
    allmargin_woComp = [];
end