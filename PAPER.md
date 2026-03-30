_Formulazione Matematica, Integratori Numerici e Visualizzazione_  
**[Claudio Simeoni / Filippo Cioccolanti]**

---
## Indice

1. Panoramica del Progetto
2. Formulazione del problema
3. Formalismo matematico
4. Scelte implementative
5. Metodi Numerici
6. Risultati della Simulazione

---
## Panoramica del Progetto

L’obiettivo di questo progetto è simulare l'evoluzione del Sistema Solare in un periodo di 20 anni. Per valutare l'entità dell'errore commesso dai vari metodi numerici utilizzati confrontiamo i nostri risultati con i dati pubblici ufficiali sul sistema solare forniti dalla NASA.

Come formalismo ci basiamo sulla gravitazione Newtoniana e sulla meccanica Hamiltoniana.

Per realizzare le simulazioni abbiamo usato il linguaggio Python, con l'ausilio delle librerie`numpy, matplotlib e astropy`.

L'intero codice del progetto è disponibile al link github https://github.com/claudioSimeoni/physics_simulations.git . Sebbene in questo documento analizziamo solo la simulazione relativa al sistema solare, il nostro codice è ben più generico e permette di trattare sistemi molto diversi, nonché visualizzare i risultati tramite animazioni. 

---
## Formulazione del problema

Il problema degli N-corpi studia il moto in tre dimensioni di un sistema formato da N corpi che interagiscono tra loro, nel nostro caso secondo il modello gravitazionale Newtoniano. 
Dato che è dimostrabilmente impossibile risolvere analiticamente questo problema mediante quadrature (dunque usando funzioni elementari) le uniche soluzioni possibili sono di carattere numerico, ragion per cui il problema è stato storicamente di grande importanza nel campo delle simulazioni numeriche.

---
## Formalismo matematico
#### Modello Hamiltoniano

Per utilizzare le equazioni di Hamilton rappresentiamo le posizioni dei pianeti mediante un vettore $\mathbf{q} \in \mathbb{R}^{3N}$ e i momenti degli stessi tramite un vettore $\mathbf{p} \in \mathbb{R}^{3N}$.
Le equazioni differenziali del sistema possono quindi essere espresse come:

$$\dot{\mathbf{q}} = \frac{\partial \mathcal{H}}{\partial \mathbf{p}}, \qquad \dot{\mathbf{p}} = - \frac{\partial \mathcal{H}}{\partial \mathbf{q}}.$$

dove $\mathcal{H}(\mathbf{q}, \mathbf{p})$ è l'hamiltoniana del sistema.

#### Derivazione dell'Hamiltoniana

L'Hamiltoniana di questo sistema è data dalla somma dell'energia cinetica $T(\mathbf{p})$ e dell'energia potenziale gravitazionale $U(\mathbf{q})$, quindi:

$$\mathcal{H}(\mathbf{q}, \mathbf{p}) = T(\mathbf{p}) + U(\mathbf{q}).$$

L'energia cinetica è data da:

$$T(\mathbf{p}) = \sum_{i=1}^N \frac{1}{2 m_i} \|\mathbf{p}_i\|^2$$

mentre l'energia potenziale gravitazionale da:


$$U_{ij} = - \frac{G m_i m_j}{\|\mathbf{q}_i - \mathbf{q}_j\|}, \quad U(\mathbf{q}) = - \sum_{1 \le i < j \le N} \frac{G m_i m_j}{\|\mathbf{q}_i - \mathbf{q}_j\|}.$$

dunque l'espressione totale dell'Hamiltoniana in termini di $\mathbf{p}$ e $\mathbf{q}$ è:


$$\mathcal{H}(\mathbf{q},\mathbf{p}) = \sum_{i=1}^N \frac{1}{2 m_i} \|\mathbf{p}_i\|^2 \quad - \sum_{1 \le i < j \le N} \frac{G m_i m_j}{\|\mathbf{q}_i - \mathbf{q}_j\|}. $$

#### Derivazione delle Equazioni di Hamilton

A questo punto è semplice scrivere le equazioni differenziali del sistema in forma chiusa:

$$\dot{\mathbf{q}_i} = \frac{\partial H}{\partial \mathbf{p}_i} = \frac{\mathbf{p}_i}{m_i}$$
$$\dot{\mathbf{p}_i} = -\frac{\partial H}{\partial \mathbf{q}_i} = G \sum_{\substack{j=1 \\ j \neq i}}^N m_i m_j \frac{\mathbf{q}_i - \mathbf{q}_j}{\|\mathbf{q}_i - \mathbf{q}_j\|^3}.$$

---
## Scelte implementative

In questa sezione discutiamo alcune scelte di carattere pratico che abbiamo effettuato nella pianificazione del codice per il progetto.

#### Forza gravitazionale e softening

Nelle simulazioni software del problema è molto comune introdurre un fattore $\varepsilon$, detto **softening**, che viene sommato alla distanza calcolata tra due corpi al fine di evitare che una collisione porti a divisioni per zero e ad errori significativi. 
Il potenziale gravitazionale modificato assume quindi questo aspetto:

$$U(\mathbf{q}) = - \sum_{1 \le i < j \le N} \frac{G m_i m_j}{\sqrt{\|\mathbf{q}_i - \mathbf{q}_j\|^2 + \varepsilon^2}}.$$

Dato che in una simulazione accurata del sistema solare che coinvolga solo i pianeti principali non ci sono collisioni abbiamo deciso di assegnare ad $\varepsilon$ un valore nullo, per ottenere risultati più vicini a quelli reali.

#### Utilizzo di Python 

Abbiamo scelto di utilizzare il linguaggio di programmazione Python, rispetto ad altre opzioni quali Matlab, valutando diversi aspetti:
- La qualità delle librerie di calcolo numerico e visualizzazione, come `numpy` e `matplotlib`;
- La struttura *object-oriented* del linguaggio che permette di scrivere codice pulito, modulare e di semplice organizzazione;
- La presenza di librerie quali `astropy`, che ci hanno permesso di interagire in modo semplice con i dataset del sistema solare mantenuti dalla NASA;

#### Utilizzo di Numpy

Numpy è la più conosciuta, affidabile e completa libreria di Python per gestire calcoli numerici complessi tra vettori e matrici. L'abbiamo scelta perché possiede una serie di proprietà che rendono le operazioni più comuni in questo ambito semplici da scrivere e veloci nell'esecuzione. Una di queste è sicuramente il *broadcasting*, che permette di effettuare operazioni aritmetiche tra array di diverse dimensioni, semplificando di molto diversi calcoli nel nostro codice.
Al fine di utilizzare al meglio le potenzialità della libreria, abbiamo optato per la rappresentazione delle coordinate del sistema dinamico:

$$\mathbf{q} = (x_1,\dots,x_N,\, y_1,\dots,y_N,\, z_1,\dots,z_N),$$
$$\mathbf{p} = (p_{x,1},\dots,p_{x,N},\, p_{y,1},\dots,p_{y,N},\, p_{z,1},\dots,p_{z,N}).$$

---
## Metodi Numerici

Abbiamo confrontato i risultati dei seguenti metodi numerici:

- **Eulero Esplicito:** metodo del primo ordine, semplice nella formulazione, non simplettico.
- **Euler Simplettico:** metodo del primo ordine simplettico.
- **Störmer-Verlet:** metodo del secondo ordine simplettico.
- **Runge-Kutta 4**: metodo del quarto ordine non simplettico, offre precisione maggiore in cambio di un calo di prestazioni.

---
## Risultati della Simulazione

Abbiamo simulato il sistema con diversi passi temporali, compresi tra 1 giorno e 0.01 giorni (~15 minuti) per un tempo totale di 10 anni. Abbiamo poi generato diversi grafici per analizzare i dati ottenuti.

### Analisi della conservazione dell'energia

In questa sezione analizziamo la conservazione dell'energia totale del sistema lungo un intervallo temporale di 20 anni, al variare del timestep e del metodo numerico utilizzato. Per ciascun grafico consideriamo l'errore relativo percentuale sull'energia, definito come $\epsilon_r(t) = \frac{H(t) - H_0}{H_0} \cdot 100$, dove $H_0$ rappresenta l'energia iniziale del sistema. In questo modo è possibile valutare quanto ciascun metodo si discosti nel tempo dalla conservazione ideale dell'energia.

---
#### Timestep = 1 giorno

Nel primo grafico sono riportati i risultati ottenuti con tutti e quattro i metodi numerici considerati.

![[docs/energy1.png]](docs/energy1.png)

Poiché l'errore dell'Eulero esplicito risulta molto maggiore rispetto agli altri metodi, nel grafico seguente lo rimuoviamo per osservare meglio il comportamento degli altri integratori.

![[docs/energy1noeuler.png]](docs/energy1noeuler.png)

Per rendere ancora più leggibile il confronto tra i metodi con errore minore, mostriamo infine un ulteriore grafico in cui viene rimosso anche l'Eulero simplettico.

![[docs/energy1noeulernosympleuler.png]](docs/energy1noeulernosympleuler.png)

Gli errori relativi massimi sull'energia, dopo una simulazione di 20 anni, sono i seguenti:

- Errore di Eulero esplicito: $4.451 \%$  
- Errore di Eulero simplettico: $0.011611 \%$  
- Errore di Verlet-Störmer: $0.000277 \%$  
- Errore di RK4: $1.648 \cdot 10^{-5} \%$

---
#### Timestep = 0.1 giorni

Per le stesse considerazioni di prima mostriamo 3 distinti grafici.

![[docs/energy01.png]](docs/energy01.png)

![[docs/energy01noeuler.png]](docs/energy01noeuler.png)

![[docs/energy01noeulernosympleuler.png]](docs/energy01noeulernosympleuler.png)

Gli errori relativi massimi sull'energia, dopo 20 anni di simulazione, risultano:

- Errore di Eulero esplicito: $1.295 \%$
- Errore di Eulero simplettico: $7.986 \cdot 10^{-4} \%$
- Errore di Verlet-Störmer: $2.620 \cdot 10^{-6} \%$
- Errore di RK4: $1.526 \cdot 10^{-10} \%$

---
#### Timestep = 0.01 giorni

Per le stesse considerazioni di prima mostriamo 3 distinti grafici.

![[docs/energy001.png]](docs/energy001.png)

![[docs/energy001noeuler.png]](docs/energy001noeuler.png)

![[docs/energy001noeulernosympleuler.png]](docs/energy001noeulernosympleuler.png)

Gli errori relativi massimi sull'energia dopo 20 anni sono:

- Errore di Eulero esplicito: $0.281 \%$  
- Errore di Eulero simplettico: $0.000106 \%$  
- Errore di Verlet-Störmer: $2.802 \cdot 10^{-8} \%$  
- Errore di RK4: $1.604 \cdot 10^{-11} \%$

---
#### Considerazioni globali

Dai risultati ottenuti emerge chiaramente che i diversi metodi numerici presentano comportamenti molto differenti dal punto di vista della conservazione dell'energia.

L'Eulero esplicito è il metodo che mostra le prestazioni peggiori: l'errore energetico cresce in modo significativo e rimane di diversi ordini di grandezza superiore rispetto agli altri metodi anche riducendo il timestep. Questo è coerente con il fatto che si tratta di un metodo del primo ordine non simplettico, che tende a introdurre una variazione considerevole dell'energia nel lungo periodo.

L'Eulero simplettico, pur essendo anch'esso un metodo del primo ordine, mostra un comportamento nettamente migliore. La ragione principale è che si tratta di un metodo simplettico: pur non conservando esattamente l'energia a ogni passo, riesce a preservare meglio la struttura hamiltoniana del sistema, mantenendo l'errore energetico limitato e tipicamente oscillante. Dai dati si osserva infatti che, al diminuire del timestep, l'errore finale decresce in modo approssimativamente lineare, in accordo con l'ordine del metodo.

Il metodo di Verlet-Störmer fornisce risultati ancora migliori. Essendo un metodo del secondo ordine e simplettico, combina una maggiore accuratezza locale con una buona conservazione della struttura geometrica del problema. Nei risultati numerici si osserva infatti una riduzione dell'errore energetico compatibile con un andamento quadratico rispetto al timestep, come previsto teoricamente.

Il metodo RK4 risulta essere quello con l'errore finale più piccolo in tutti i casi considerati. Tuttavia, è importante sottolineare che RK4, pur essendo un metodo di ordine quattro, non è simplettico. Di conseguenza, l'ottima accuratezza osservata sull'energia finale non implica necessariamente una conservazione strutturale dell'energia nel lungo periodo, ma riflette soprattutto l'elevata precisione del metodo per timestep sufficientemente piccoli. Inoltre, per valori molto piccoli di $h$, l'errore finale sull'energia diventa così ridotto da poter essere influenzato anche da effetti di arrotondamento numerico e cancellazione, motivo per cui non si osserva sempre una legge di scala perfettamente riconducibile all'ordine teorico del metodo.

### Analisi delle traiettorie dei pianeti

Nei seguenti grafici mettiamo a confronto la performance dei 4 metodi numerici nel calcolare le orbite di Mercurio, Terra, Giove e Nettuno. Abbiamo scelto di presentare questi pianeti per analizzare sia pianeti del sistema solare interno che esterno, con orbite più o meno eccentriche.
Nei grafici viene mostrato nell'asse delle ordinate l'errore assoluto in km tra la posizione effettiva del pianeta considerato e quella da noi calcolata.

##### Primo sguardo ai pianeti

![[docs/earth1.png]](docs/earth1.png) *Terra* 

![[docs/mercury1.png]](docs/mercury1.png)
*Mercurio*

![[docs/jupiter1.png]](docs/jupiter1.png)
*Giove*

![[docs/neptune1.png]](docs/neptune1.png)
*Nettuno*

In questi primi mettiamo a confronto i quattro pianeti con timestep di 1 giorno. Notiamo subito che Eulero esplicito ha prestazioni molto peggiori degli altri metodi, specialmente nel caso dei pianeti del sistema solare interno. In questo metodo, gli errori allontanano sempre i pianeti dalle loro orbite e dunque l'accumularsi di queste imprecisioni tutte dello stesso tipo causa un errore esponenziale nel tempo. 
I metodi simplettici, d'altro canto, preservano i volumi nello spazio delle fasi e dunque costringono le orbite ad oscillare di poco intorno a quelle vere, garantendo stabilità a lungo termine. Può dunque risultare insolito che RK4 commetta un errore così piccolo nonostante non sia simplettico, ma ciò è dovuto alla sua natura di metodo di quarto ordine: su scale di tempo così brevi non è possibile osservare un errore sensibile per la grande precisione locale di questo metodo, mentre per Eulero è più che sufficiente. 

##### Differenze tra pianeti interni ed esterni

Interessante è anche la differenza di comportamento degli integratori tra i pianeti del sistema solare interno ed esterno. Guardando ad esempio il grafico di Nettuno, la differenza di prestazione tra Eulero e gli altri metodi diminuisce, così come l'entità dell'errore (minore di circa 2 ordini di grandezza). Questo è dovuto a due fattori principali:
- la distanza di Nettuno dal Sole, che rende perturbazioni sulla sua posizione quasi ininfluenti al fine della corretta determinazione della forza gravitazionale su di esso. Per Mercurio o per la Terra, invece, piccoli errori sono molto più significativi;
- la minor velocità di Nettuno, che rende minore l'errore di posizione locale e conseguentemente l'errore nel calcolo delle forze.

##### Analisi globale dei dati raccolti

Di seguito le tabelle con tutti i dati raccolti, divisi per timestep, che ci saranno comode per fare alcune osservazioni.

| **Pianeta**  | **Eulero Esplicito** | **Eulero Simplettico** | **Verlet Stormer** | **RK4**  |
| ------------ | -------------------- | ---------------------- | ------------------ | -------- |
| **Mercurio** | 69.805.515,61        | 1.533.382,29           | 5.618,39           | 721,96   |
| **Terra**    | 289.848.671,46       | 80.211,13              | 930,85             | 762,53   |
| **Giove**    | 1.874.086,00         | 10.351,65              | 3.224,86           | 3.224,31 |
| **Nettuno**  | 2.862,19             | 1.605,62               | 541,09             | 541,09   |
*Timestep 0.01 giorni* 

| **Pianeta**  | **Eulero Esplicito** | **Eulero Simplettico** | **Verlet Stormer** | **RK4**  |
| ------------ | -------------------- | ---------------------- | ------------------ | -------- |
| **Mercurio** | 203.498.521,20       | 14.430.640,84          | 511.826,73         | 707,04   |
| **Terra**    | 118.513.926,05       | 766.594,69             | 19.531,11          | 761,44   |
| **Giove**    | 18.667.044,82        | 114.521,13             | 3.279,40           | 3.224,31 |
| **Nettuno**  | 23.801,33            | 20.319,63              | 541,11             | 541,09   |
*Timestep 0.1 giorni*

| **Pianeta**  | **Eulero Esplicito** | **Eulero Simplettico** | **Verlet Stormer** | **RK4**      |
| ------------ | -------------------- | ---------------------- | ------------------ | ------------ |
| **Mercurio** | 690.694.221,29       | 78.903.493,93          | 49.724.647,24      | 1.671.018,46 |
| **Terra**    | 505.797.767,03       | 4.192.980,76           | 1.864.538,57       | 8.317,65     |
| **Giove**    | 181.724.848,82       | 1.150.239,16           | 9.147,57           | 3.224,31     |
| **Saturno**  | 20.309.121,95        | 1.580.231,14           | 3.841,89           | 2.867,18     |
| **Nettuno**  | 233.203,19           | 207.607,35             | 543,42             | 541,09       |
*Timestep 1 giorno*

La principale osservazione attuabile osservando queste tabelle è l'aderenza dei metodi al loro ordine teorico. Confrontando ad esempio le prestazioni di RK4 nel calcolo dell'orbita di Mercurio con timestep di 1 giorni e 0.1 giorni, ci rendiamo conto che il rapporto tra i due errori è di circa $10^4$ , coerentemente con la natura di RK4 di metodo di quarto ordine.
Similmente, se compariamo le prestazioni di Eulero Simplettico e Stormer-Verlet sull'orbita di Nettuno con timestep di 0.1 giorni notiamo che il rapporto è di circa $10^2$, coerentemente con le loro nature di metodo di primo e secondo ordine rispettivamente.
Questa coerenza con la teoria ha comunque un limite, probabilmente dato da problemi di precisione relativi all'aritmetica con virgola mobile: dalle tabelle è ben chiaro che una volta raggiunta una certa precisione anche diminuendo il timestep non si ha un miglioramento significativo delle precisioni. A questo proposito si osservi Nettuno, che quando simulato con RK4 o Stormer-Verlet mostra sempre lo stesso errore di posizione.

---

## Conclusioni

Dal nostro studio emerge la netta superiorità di Runge-Kutta 4 rispetto agli altri metodi, nonostante la sua assenza di simpletticità. La scala temporale considerata è infatti troppo breve rispetto alla precisione del metodo perché questa proprietà diventi influente sulla sua precisione. Nelle simulazioni con timestep 0.01 giorni, il più piccolo tra quelli considerati, l'entità dell'errore su Giove e Nettuno è veramente piccola, dell'ordine del ~3% del loro raggio, mentre sulla Terra e su Mercurio è maggiore data la maggiore velocità della loro orbita. Ci riteniamo soddisfatti dei risultati ottenuti, dati i mezzi e la portata del progetto.