function bonobo=GridIndexFinding(bonobo,Grid)

    nObj=numel(bonobo.Cost);
    
    nGrid=numel(Grid(1).LB);
    
    bonobo.GridSubIndex=zeros(1,nObj);
    
    for j=1:nObj
        
        bonobo.GridSubIndex(j)=...
            find(bonobo.Cost(j)<Grid(j).UB,1,'first');
        
    end

    bonobo.GridIndex=bonobo.GridSubIndex(1);
    for j=2:nObj
        bonobo.GridIndex=bonobo.GridIndex-1;
        bonobo.GridIndex=nGrid*bonobo.GridIndex;
        bonobo.GridIndex=bonobo.GridIndex+bonobo.GridSubIndex(j);
    end
    
end