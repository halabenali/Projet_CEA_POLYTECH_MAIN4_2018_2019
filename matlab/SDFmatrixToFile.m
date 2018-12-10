function SDFmatrixToFile(file_name, seq_starter, seq_increment, seq_end)

% Ecrit dans le fichier file_name des matrices symetriques definies
% positives de dimension (n,n). Pour n allant de seq_starter a seq_end par
% pas de seq_increment (n = seq_starter : seq_increment : seq_end)
% Le fichier file_name resultant est place dans le repertoire
% ../data/matrices

% ATTENTION : Si le fichier existe deja dans le repertoire ../data/matrices
% son contenu est ecrase et remplace par les matrices generees

% Tous les arguments doivent etre des chaines de caracteres ou string

file_name = "../data/matrices/"+file_name;

%% creation/ouverture du fichier et ecriture de l'en tete de la premiere matrice
fid = fopen(file_name,'wt');
fprintf(fid, '%s\n', "M " + seq_starter);
fclose(fid);

%% conversion des arguments numeriques
seq_starter = str2num(seq_starter);
seq_increment = str2num(seq_increment);
seq_end = str2num(seq_end);

%% ecriture de la premiere matrice dans le fichier

% generation de la matrice
val_propres = rand(seq_starter,1);
M = sprandsym(seq_starter,rand(),val_propres);
M = full(M);
A = M;

% ecriture
dlmwrite(file_name,M,'-append','delimiter','\t','precision',3);

%% ecriture des autres matrices 
for i = seq_starter + seq_increment : seq_increment : seq_end
    % ecriture de l'en tete
    fid = fopen(file_name,'at');
    fprintf(fid, '%s\n', "M " + num2str(i,'%d'));
    fclose(fid);
    % generation de la matrice
    val_propres = rand(i,1);
    M = sprandsym(i,rand(),val_propres);
    M = full(M);
    % ecriture de la matrice
    dlmwrite(file_name,M,'-append','delimiter','\t','precision',3);
end

end