time = (0:h)';    % time horizon for IRFs and FEVDs

figure('Position',[100 100 1000 600],'PaperPositionMode','Auto','DefaultAxesFontSize',13);
for j=1:nvar 
    h(j) = subplot(ceil(nvar/3),3,j);

    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper_pe(1,j); IRFslower_pe(1:end,j); flipud([IRFsupper_pe(1:end,j); IRFslower_pe(end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.2);
    set(hh,'edgecolor','none'); 

    hold on;
    hh=fill([time(1); time(1:end); flipud([time(1:end); time(end)])],[IRFsupper2_pe(1,j); IRFslower2_pe(1:end,j); flipud([IRFsupper2_pe(1:end,j); IRFslower2_pe(end,j)])],[0.1, 0.4470, 0.7410]); 
    set(hh,'facealpha',.4);
    set(hh,'edgecolor','none');

    p1=plot(time, IRFs_pe(:,j),'k', 'Linewidth', 1.5); 
    if ~ismember(0,get(gca,'ylim'))
        line(get(gca,'xlim'),[0 0],'Color','k')
    end
    box on
    grid on ;hold off;
    title(varNames_paper{j}) 
    ylabel('\%');
    xlim([0,horizon]);
    if dataFrequency == 'M'
        xlabel('Months');
        xticks(0:10:horizon)
    elseif dataFrequency == 'Q'
        xlabel('Quarters');
        xticks(0:5:horizon)
    end   
end
pause(0.001)
h=axes('Position',[0.25,0,.5,.5],'Xlim',[0 1],'Ylim',[0 1]);
set(h,'Visible','off');
if k==1 && strcmp(estType,'proxy')
    string_1stage = ['First stage regression: F: ',num2str(olsEst.r1.F,' %2.2f'),', robust F: ',num2str(olsEst.r1.Frobust,' %2.2f'),', $R^2$: ',num2str(olsEst.r1.R2*100,' %1.2f'),'\%, Adjusted $R^2$: ',num2str(olsEst.r1.R2adj*100,' %1.2f'),'\%'];
    text('Position',[-0.16 -0.002],'string',string_1stage,'FontSize',14);
end
tightfig;
if saveFigs
    print('-dpdf', gcf, strcat('figures/IRFs_',estType));  
end