function rep=RepMemberDeletion(rep,gamma)

    % Grid Index of All Repository Members
    GI=[rep.GridIndex];
    
    % Occupied Cells
    OC=unique(GI);
    
    % Number of Bonobos in Occupied Cells
    N=zeros(size(OC));
    for k=1:numel(OC)
        N(k)=numel(find(GI==OC(k)));
    end
    
    % Selection Probabilities
    P=exp(gamma*N);
    P=P/sum(P);
    
    % Selected Cell Index using Roulette Wheel Selection method
    sci=RWSelection(P);
    
    % Selected Cell
    sc=OC(sci);
    
    % Selected Cell Members
    SCM=find(GI==sc);
    
    % Selected Member Index
    smi=randi([1 numel(SCM)]);
    
    % Selected Member
    sm=SCM(smi);
    
    % Delete Selected Member
    rep(sm)=[];

end