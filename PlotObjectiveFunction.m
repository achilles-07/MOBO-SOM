function PlotObjectiveFunction(pop,rep,nObj)
if(nObj==2)
    pop_costs=[pop.Cost];
    plot(pop_costs(1,:),pop_costs(2,:),'ko');
    hold on;
    rep_costs=[rep.Cost];
    plot(rep_costs(1,:),rep_costs(2,:),'r*');
    xlabel('1^{st} Objective');
    ylabel('2^{nd} Objective');
    grid on;
    hold off;
   
elseif (nObj==3)
    pop_costs=[pop.Cost];
    plot3(pop_costs(1,:),pop_costs(2,:),pop_costs(3,:), 'ko');
    hold on;
    rep_costs=[rep.Cost];
    plot3(rep_costs(1,:),rep_costs(2,:),rep_costs(3,:),'r*');
    xlabel('1^{st} Objective');
    ylabel('2^{nd} Objective');
    zlabel('3^{rd} Objective');
    grid on;
    hold off;
else
    disp('2D or 3D Graph not possible');
end

end