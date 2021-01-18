## genepop-shiny

Le projet genepop-shiny utilise shiny afin d'offrir une interface au package [genepop](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop)

## Architecture

```
├── genepop-shiny/
│   ├── app.R
│   ├── R/
│   │   ├── menugauche.R
│   |   ├── helper_functions.R  
│   ├── pages/
│   ├── opts/
│   ├── www/
│   ├── README.md
│   ├── config.yml
│   ├── Dockerfile
│   ├── application.yml
```

Le fichier app.R définit la structure globale de l'apppli, avec ses deux composantes :
* UI pour l'aspect affichage des différentes pages
* server pour l'aspect réaction aux différentes actions sur les composants de l'interface.

Le design fait appel à shinydashboard.

Les menus sont définis dans le fichier menugauche.R.

Chaque menu principal est défini par une fonction menuItem et chaque sous menu par une fonction menuSubItem.

Pour chaque menu on défini le nom de la page (tabName) à charger suite à un clic.

La définition des composants des pages se fait dans les fichiers pages/pages_def_* (regroupés par option de genepop)

La définion de la logique des réactions lors des modification des composants de chaque page est dans les fichiers opts/opt**.R

## Dépendance R
* [shiny](https://shiny.rstudio.com/)
* [shinydashboard](https://rstudio.github.io/shinydashboard/)
* [shinyjs](https://cran.r-project.org/web/packages/shinyjs/index.html)
* [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
* [genepop](http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop)
* [yaml](https://cran.r-project.org/web/packages/yaml/index.html)

## Personnalisation

Le fichier config.yml est à votre disposition afin de personnaliser genepop-shiny.<br/>

## Utilisation

* Première étape : Dans l'onglet "Load Data", vous devez charger votre fichier de données au format genepop. Vous pouvez le voir dans le style brut ou dans une table, dans la partit "File content".

* Deuxième étape : Rendez-vous dans une des options afin de réaliser une des fonctionnalités de genepop. Après une sélection des paramètres d'entrée, utilisez le bouton "Run" pour lancer les calculs.

* Troisième étape : une fois les calculs terminés vous pouvez visualiser les résultats, ainsi que les télécharger à l'aide du bouton "Dowload".
