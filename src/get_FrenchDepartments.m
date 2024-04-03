% List all the French departments, including overseas territories in the 
% order they are stored in the SIDEP repository available on
% https://www.data.gouv.fr/fr/datasets/donnees-de-laboratoires-pour-le-depistage-a-compter-du-18-05-2022-si-dep/.
%
% B. Pascal,
% March, 2024


function AllDepartments = get_FrenchDepartments(display)

    % Inputs:  - display: 1 for displaying the list of all available countries,
    %                     0 does noting.
    %
    % Ouputs:  - AllDepartments: structure containing:
    %                     Names: list of French departments names
    %                     Number: official INSEE numerotation of French departments
    %                     Indices: index of each department in the storage from 1 to D
    %
    % Nota bene: Corsica is divided into two departments, with numerotation
    % 2A and 2B and stored in positions 30 and 31.

    % Number of French departments including overseas territories
    D           = 104;

    % Names of the French departments for displaying the table
    Names_Table = ["Ain",
        "Aisne",
        "Allier",
        "Alpes-de-Haute-Provence",
        "Hautes-Alpes",
        "Alpes-Maritimes",
        "Ardèche",
        "Ardennes",
        "Ariège",
        "Aube",
        "Aude",
        "Aveyron",
        "Bouches-du-Rhône",
        "Calvados",
        "Cantal",
        "Charente",
        "Charente-Maritime",
        "Cher",
        "Corrèze",
        "Côte-d'Or",
        "Côtes-d'Armor",
        "Creuse",
        "Dordogne",
        "Doubs",
        "Drôme",
        "Eure",
        "Eure-et-Loir",
        "Finistère",
        "Corse-du-Sud",
        "Haute-Corse",
        "Gard",
        "Haute-Garonne",
        "Gers",
        "Gironde",
        "Hérault",
        "Ille-et-Vilaine",
        "Indre",
        "Indre-et-Loire",
        "Isère",
        "Jura",
        "Landes",
        "Loir-et-Cher",
        "Loire",
        "Haute-Loire",
        "Loire-Atlantique",
        "Loiret",
        "Lot",
        "Lot-et-Garonne",
        "Lozère",
        "Maine-et-Loire",
        "Manche",
        "Marne",
        "Haute-Marne",
        "Mayenne",
        "Meurthe-et-Moselle",
        "Meuse",
        "Morbihan",
        "Moselle",
        "Nièvre",
        "Nord",
        "Oise",
        "Orne",
        "Pas-de-Calais",
        "Puy-de-Dôme",
        "Pyrénées-Atlantiques",
        "Hautes-Pyrénées",
        "Pyrénées-Orientales",
        "Bas-Rhin",
        "Haut-Rhin",
        "Rhône",
        "Haute-Saône",
        "Saône-et-Loire",
        "Sarthe",
        "Savoie",
        "Haute-Savoie",
        "Paris",
        "Seine-Maritime",
        "Seine-et-Marne",
        "Yvelines",
        "Deux-Sèvres",
        "Somme",
        "Tarn",
        "Tarn-et-Garonne",
        "Var",
        "Vaucluse",
        "Vendée",
        "Vienne",
        "Haute-Vienne",
        "Vosges",
        "Yonne",
        "Territoire de Belfort",
        "Essonne",
        "Hauts-de-Seine",
        "Seine-Saint-Denis",
        "Val-de-Marne",
        "Val-d'Oise",
        "Guadeloupe",
        "Martinique",
        "Guyane",
        "La Réunion",
        "Saint-Pierre-et-Miquelon",
        "Mayotte",
        "Saint-Barthélemy",
        "Saint-Martin"];

    % Names of the French departments for displaying the table
    Names = ["Ain",
        "Aisne",
        "Allier",
        "Alpes-de-Haute-Provence",
        "Hautes-Alpes",
        "Alpes-Maritimes",
        "Ard{\`e}che",
        "Ardennes",
        "Ari{\`e}ge",
        "Aube",
        "Aude",
        "Aveyron",
        "Bouches-du-Rh{\^o}ne",
        "Calvados",
        "Cantal",
        "Charente",
        "Charente-Maritime",
        "Cher",
        "Corr{\`e}ze",
        "C{\^o}te-d'Or",
        "C{\^o}tes-d'Armor",
        "Creuse",
        "Dordogne",
        "Doubs",
        "Dr{\^o}me",
        "Eure",
        "Eure-et-Loir",
        "Finist{\`e}re",
        "Corse-du-Sud",
        "Haute-Corse",
        "Gard",
        "Haute-Garonne",
        "Gers",
        "Gironde",
        "H{\'e}rault",
        "Ille-et-Vilaine",
        "Indre",
        "Indre-et-Loire",
        "Is{\`e}re",
        "Jura",
        "Landes",
        "Loir-et-Cher",
        "Loire",
        "Haute-Loire",
        "Loire-Atlantique",
        "Loiret",
        "Lot",
        "Lot-et-Garonne",
        "Loz{\`e}re",
        "Maine-et-Loire",
        "Manche",
        "Marne",
        "Haute-Marne",
        "Mayenne",
        "Meurthe-et-Moselle",
        "Meuse",
        "Morbihan",
        "Moselle",
        "Ni{\`e}vre",
        "Nord",
        "Oise",
        "Orne",
        "Pas-de-Calais",
        "Puy-de-D{\^o}me",
        "Pyr{\'e}n{\'e}es-Atlantiques",
        "Hautes-Pyr{\'e}n{\'e}es",
        "Pyr{\'e}n{\'e}es-Orientales",
        "Bas-Rhin",
        "Haut-Rhin",
        "Rh{\^o}ne",
        "Haute-Sa{\^o}ne",
        "Sa{\^o}ne-et-Loire",
        "Sarthe",
        "Savoie",
        "Haute-Savoie",
        "Paris",
        "Seine-Maritime",
        "Seine-et-Marne",
        "Yvelines",
        "Deux-S{\`e}vres",
        "Somme",
        "Tarn",
        "Tarn-et-Garonne",
        "Var",
        "Vaucluse",
        "Vend{\'e}e",
        "Vienne",
        "Haute-Vienne",
        "Vosges",
        "Yonne",
        "Territoire de Belfort",
        "Essonne",
        "Hauts-de-Seine",
        "Seine-Saint-Denis",
        "Val-de-Marne",
        "Val-d'Oise"
        "Guadeloupe",
        "Martinique",
        "Guyane",
        "La R{\'e}union",
        "Saint-Pierre-et-Miquelon",
        "Mayotte",
        "Saint-Barth{\'e}lemy",
        "Saint-Martin"];

    % Numerotation of French departments
    Numbers = strings(D,1);

    % Indices of the departments
    Indices = zeros(D,1);

    FrenchDepartments = cell(length(Names_Table),3);
    for d = 1:length(Names_Table)
        FrenchDepartments{d,1}     = Names_Table(d);
        if d < 20
            Numbers(d)             = sprintf("%02d",d);
            Indices(d)             = d;
            FrenchDepartments{d,2} = Numbers(d);
            FrenchDepartments{d,3} = Indices(d);
        elseif d >= 20 && d < 29
            Numbers(d)             = sprintf("%02d",d+1);
            Indices(d)             = d;
            FrenchDepartments{d,2} = Numbers(d);
            FrenchDepartments{d,3} = Indices(d);
        elseif d == 29
            Numbers(d)             = "2A";
            Indices(d)             = d;
            FrenchDepartments{d,2} = Numbers(d);
            FrenchDepartments{d,3} = Indices(d);
        elseif d == 30
            Numbers(d)             = "2B";
            Indices(d)             = d;
            FrenchDepartments{d,2} = Numbers(d);
            FrenchDepartments{d,3} = Indices(d);
        elseif d > 30 && d <= 96
            Numbers(d)             = sprintf("%02d",d-1);
            Indices(d)             = d;
            FrenchDepartments{d,2} = Numbers(d);
            FrenchDepartments{d,3} = Indices(d);
        else
            Numbers(d)             = sprintf("%03d",970 + (d - 96));
            Indices(d)             = d;
            FrenchDepartments{d,2} = Numbers(d);
            FrenchDepartments{d,3} = Indices(d);
        end
    end

    % Store names, numerotation and indices for output
    AllDepartments.Names_Table     = Names_Table; 
    AllDepartments.Names           = Names; 
    AllDepartments.Numbers         = Numbers; 
    AllDepartments.Indices         = Indices; 

    % Display the recapitulative table if required
    if display
        disp(cell2table(FrenchDepartments,'VariableNames',["Name","Number","Index"]))
    end

end