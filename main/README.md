# Simulazione

### simulazione.conf : file di configurazione
**TIMER:** indica il tempo di simulazione (numero di iterazioni fatte).

**N_VERTEX:** numero di vertici del grafo da generare.

**N_RAN:** se espresso, verranno fatte N_RAN simulazioni ed i risultati saranno mediati (per quanto riguarda le frequenze).

**WITH_STEP:** indica se applicare un'euristica diversa per la scelta del vicino in cui andare ad inserire la pallina

- **WITH_STEP = 1**, si seleziona un nodo e si scelgono **K** vicini(*entrambe le selezioni sono fatte in modo uniforme*): tra questi viene scelto il nodo che ha il minor numero di palline, eventualmente scegliendo anche se stesso (nodo selezionato all'inizio).
- **WITH_STEP = 0**, si seleziona un nodo e successivamente un vicino al quale verrà inserita la palline (*le selezioni vengono fatte in modo uniforme*).

**K**: indica il numero di nodi adiacenti (vicini) da selezionare nel caso **WITH_STEP = 1**.

**DEATH_PROB**: se indicata, essa è la probabilità con la quale una pallina muore. 

**TYPE_OF_GRAPH**: indica il tipo di grafo della simulazione:

- **FULL_GRAPH**: grafo totalmente connesso con numero di archi pari a N_VERTEX * (N_VERTEX - 1)

- **DEGREE_SEQUENCE**: grafo generato tramite sequenza di gradi. In questo caso occorrerà indicare altri due valori, **DEGREE1** e **DEGREE2** che indicano, rispettivamente, che la metà dei nodi nel grafo avrà grado *DEGREE1* e l'altra metà *DEGREE2*.

- **REGULAR_GRAPH**: grafo regolare di grado **DEGREE1** da indicare. 
NB: si può ottenere un *REGULAR_GRAPH* a partire da un *DEGREE_SEQUENCE* scegliendo uguali *DEGREE1* e *DEGREE2*