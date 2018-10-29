package MST;// Java program for Kruskal's algorithm to find Minimum
// Spanning Tree of a given connected, undirected and
// weighted graph

import java.util.*;
import java.lang.*;

class MST {

    // Total number of vertices
    int N;
    // V-> no. of vertices & E->no.of edges
    int E;
    // collection of all edges
    Edge[] edges;

    /** Constructor for MST ------------------------------------------------------------------
     */
    MST(int n, int e) {
        N = n;  // total number of vertices
        E = e; // all edges
        edges = new Edge[E];
        for (int i = 0; i < e; i++) {
//            System.out.println("current edge index in MST is: "+i);
            edges[i] = new Edge(0,0,0);
        }
    }

    /** Inner Class for Vertices ------------------------------------------------------------------
     */
    private static class Vertices{
        ArrayList<Position> pos = new ArrayList<>();
        class Position{
            float x;
            float y;
            int n;
            Position(int n, float x, float y){
                this.n = n;
                this.x = x;
                this.y = y;
            }
        }
        void addPos(int n){
            Random rand = new Random();
            float x = rand.nextFloat();
            float y = rand.nextFloat();
            Position p = new Position(n, x, y);
            pos.add(p);
        }

        public float getPosX(int index){
            return pos.get(index).x;
        }
        public float getPosY(int index){
            return pos.get(index).y;
        }
    }

    /** Inner Class for Edge ----------------------------------------------------------------------
     */
    static class Edge implements Comparable<Edge> {
        int src, dest;
        double weight;

        public Edge(int src, int dest, double weight){
            this.dest = dest;
            this.src = src;
            this.weight = weight;
        }

        @Override
        public int compareTo(Edge b)
        {
            return this.weight > b.weight ? 1 : -1;
        }

        @Override
        public String toString(){
            return "( from  " + src  + " to " + dest + ", weight " + weight +" )";
        }

    };

    /** Inner Class for Union_Find ----------------------------------------------------------------
     */
    static class Disjoint_Set{
        private int[] root, rank;
        int num;

        Disjoint_Set(int num){
            rank = new int[num];
            root = new int[num];
            this.num = num;
            // initialization for this data structure
            init_set(num);
        }

        void init_set(int num){
            for (int i = 0; i<num; i++){
                // initialize the root to itself
                root[i] = i;
                rank[i] = 0;
            }
        }

        /** Path Compression technique ----------------------------------------------------------------
         */
        public int find_set(int a){
            if(root[a]!=a){
                root[a] = find_set(root[a]);
            }
            return root[a];
        }

        /** Union by Rank -----------------------------------------------------------------------------
         */
        public void Union(int x,int y){
            int xroot = find_set(x);
            int yroot = find_set(y);
            if(xroot == yroot){
                return;
            }
            if (rank[xroot] < rank[yroot])
                root[xroot] = yroot;

            else if(rank[xroot] > rank[yroot])
                root[yroot] = xroot;

            else{
                root[yroot] = xroot;
                rank[xroot]++;
            }
        }
    }

    /** Generate MST by Kruskal's algorithm, return the average weight -----------------------------------------------------------------------------
     */
    double KruskalMST()
    {
        // This will store the resulting MST
        Edge[] mst= new Edge[N];
//        System.out.println("the value of N is: "+ N);

        for (int i = 0; i< N-1; i++)
            // to store the final mst
            mst[i] = new Edge(0,0,0);

        // sorting the edges, num == N(N-1)/2 for complete graph
//        System.out.println("length of edges:  "+ edges.length);
        Arrays.sort(edges);

        // Object for Disjoint_Set ------------------------------------------------------
        Disjoint_Set sets = new Disjoint_Set(N);

        int i = 0;
        int e = 0;

        while (i < N-1) {


//            System.out.println(edges[e].toString());

            Edge current_edge = edges[e];

            int xroot = sets.find_set(current_edge.src);
            int yroot = sets.find_set(current_edge.dest);

                // Cycle Detection ----------------------------------------------------------
            if (xroot != yroot) {
                // not cycle -> pick
                mst[i] = current_edge;
                sets.Union(xroot, yroot);
                i++;
            }
            e++;
            }


        double total_W = 0;
        for (int m = 0; m < i; m++) {
            total_W += mst[m].weight;}
        return total_W;
    }

    PriorityQueue<Edge> Q =new PriorityQueue();
    // Construct MST by Prim's algorithm, return the average weight

    double PrimMST(){
        Edge[] Pmst= new Edge[N];
        int e = 0;
        while( !Q.isEmpty()){

//            u = Q.poll();
//
//            for(Vertex v: Vertexes){
//                if (v.label != u.label){
//                    if(Q.Contain(v) && matrix[u.label][v.label]< v.key){
//                        v.p = u;
//                        v.key = matrix[u.label][v.label];
//                    }
//                }
//            }
//
//            Q.Update();
////                System.out.println("");
//            if (u.p != null) {
//                A[index] = new Edge(new Point(u.p.label), new Point(u.label), u.key);
////                    System.out.println(A[index]);
//                index++;
//            }

        }

        double total_W = 0;

        for (Edge t : Pmst){
            total_W += t.weight;
        }
        return total_W;
    }


    /** Main Program --------------------------------------------------------------------------
     */
    public static void main (String[] args)
    {
        double[] averageWeight_Q1a = new double[5];
        long[] averageRunTime_Q1a = new long[5];

        double[] K_averWeight_Q1c = new double[5];
        long[] K_averRunTime_Q1c = new long[5];

        double[] P_averWeight_Q1c = new double[5];
        long[] P_averRunTime_Q1c = new long[5];

        int sizes[] = {100,500,1000};

//        /** Q1_a Implementation begins----------------------------------------------------------
//         */
//        for(int s =0; s<sizes.length;s++){
//
//            double testWeight_Q1a = 0.0;
//            long testTime_Q1a = (long) 0.0;
//
//
//            for (int d=0; d<50;d++) {
//                //System.out.println(d);
//                int N = sizes[s]; // Input number of vertices
//                int E = N * (N - 1) / 2; // Number of edges in complete graph
//                MST graph = new MST(N, E);
//
//                // Initialize the ArrayList of Points Objects
//                Vertices V = new Vertices();
//
//                // Assign the position for all vertices
//                for (int k = 0; k < N; k++) {
//                    V.addPos(k);}
//
//                int t = 0; // t_th edge
//                double[][] w_matrix = new double[N][N]; // matrix to store
//
//                // Q1_a Complete Graph Initialization
//                // Fill the matrix to represent the graph
//                for (int i = 0; i < N; i++) { // src = i
//                    for (int j = i; j < N; j++) { // dest = j
//                        if (i == j) {
//                            w_matrix[i][j] = 999;}
//                        else {
//                            graph.edges[t].src = i;
//                            graph.edges[t].dest = j;
//                            double w = Math.sqrt(Math.pow(V.getPosX(i) - V.getPosX(j), 2) + Math.pow(V.getPosY(i) - V.getPosY(j), 2));
//                            graph.edges[t].weight = w;
//                            w_matrix[i][j] = w; // w_matrix[src][dest] = weight
//                            t++;}
//                    }
//                }
//
//                // Timing for Q1_a
//                long startTime_Q1a = System.nanoTime();
//                testWeight_Q1a+=graph.KruskalMST();
//                long endTime_Q1a = System.nanoTime();
//                testTime_Q1a +=(endTime_Q1a - startTime_Q1a);
//
//            }
//
//            averageWeight_Q1a[s] = testWeight_Q1a/50;
//            averageRunTime_Q1a[s] =  testTime_Q1a/50;
//        }
//
//        for(int s =0; s<sizes.length;s++){
//            System.out.println("Complete Graph Size: "+ sizes[s] + " with an average Weight: "+averageWeight_Q1a[s]);
//        }
//        for(int s =0; s<sizes.length;s++){
//            System.out.println("Complete Graph Size: "+ sizes[s] + " with an average RunTime: "+averageRunTime_Q1a[s]);
//        }

        /** Q1_c Implementation begins--------------------------------------------------------------------
        */


//        Boolean unconnected = true;
        for(int s = 0; s < sizes.length; s++){

            double K_testWeight_Q1c = 0.0;
            long K_testTime_Q1c = (long) 0.0;

            double P_testWeight_Q1c = 0.0;
            long P_testTime_Q1c = (long) 0.0;

            for (int d=0; d < 50; d++) { // repeated exp

                int N = sizes[s];
                //int E = N * (N - 1) / 2; // Upper boundary of the number of edges

//                MST graph_c = new MST(N, E);

                // Initialize the ArrayList of Points Objects
                Vertices V = new Vertices();

                // Assign the position for all vertices
                for (int k = 0; k < N; k++) {
                    V.addPos(k);
                }

                /** Q1_c Random Connected Graph Initialization-------------------------------------------------
                 */
                // Adjacency list representation of graph

                LinkedList<Integer>[] AdjList = new LinkedList[N];
                ArrayList<Edge> pickedEdges = new ArrayList<Edge>();

                for (int i = 0; i < N; i++) {
                    AdjList[i] = new LinkedList<>();
                }

                // Initialize flag boolean
                Boolean unconnected = true;

                int count = 0; // count variable to track current available edges (no dup, min -> max)
                Disjoint_Set sets = new Disjoint_Set(N);

                while (unconnected) {
                    Random rand = new Random();
                    int node1_index = rand.nextInt(N);
                    int k = rand.nextInt(N);
                    while (k == node1_index) {
                        k = rand.nextInt(N); // random form  0 - (N-1)
                    }
                    int node2_index = k;

                    /** check the unique membership node2 in AdjList[node1] ------------------------------------
                     */
                    boolean unique = true;
                    if (AdjList[node1_index].contains(node2_index)){
                        unique = false;
                    }

                    if (unique) { // if this edge didn't cover, new edge to process
                        AdjList[node1_index].add(node2_index);
                        AdjList[node2_index].add(node1_index);


                        sets.Union(node1_index, node2_index);
//                        System.out.println("below is connected 2 nodes: ");
//                        System.out.println("node1_index: "+node1_index);
//                        System.out.println("node2_index: "+node2_index);

                        int min = Math.min(node1_index, node2_index);
                        int max = Math.max(node1_index, node2_index);
                        double w = Math.sqrt(Math.pow(V.getPosX(min) - V.getPosX(max), 2) + Math.pow(V.getPosY(min) - V.getPosY(max), 2));

                        Edge e = new Edge(min,max,w);
                        pickedEdges.add(e);

                        int first_root = sets.find_set(0);
                        int unconnected_src = 0;

                        /** traversal for the graph_C ---------  check whether this graph is connected
                         */
                        for (int m = 1; m < AdjList.length; m++) {
//                            System.out.println("parent of  " + m + " is "+ sets.find_set(m));
                            // if there is one group root is not the first root --> unconnected --> repeated while loop
                            if (sets.find_set(m) != first_root) {
//                                System.out.println("value of m is "+ m);

//                                System.out.println("parent of  "+ m +" is " + sets.find_set(m) + "  the first parent: "+ first_root);
//                                //System.out.println(sets.find_set(1)+ " ssssssssssss "+ sets.find_set(4));
//                                System.out.println(AdjList[0]);
//                                System.out.println(AdjList[1]);
//                                System.out.println(AdjList[2]);
//                                System.out.println(AdjList[3]);
//                                System.out.println(AdjList[4]);
                                //System.out.println("value of m is "+ m);
                                unconnected_src++;
                                // continue;
                                // the graph is still unconnected
                            }
                        }

//                            System.out.println("after checked: "+unconnected);
                        HashSet set = new HashSet();

                        for (int num = 0; num < AdjList.length; num++){
                            for (int t : AdjList[num]){
                                set.add(t);
                            }
                        }
                        int numVisited = set.size();

                        if (unconnected_src == 0 && numVisited == N) {
//                            System.out.println("hahahhahahah connected!!");
//                            // CHECK THE PARENT FOR ALL
//                            System.out.println("PARENT of 0: "+ sets.find_set(0));
//                            System.out.println("PARENT of 1: "+ sets.find_set(1));
//                            System.out.println("PARENT of 2: "+ sets.find_set(2));
//                            System.out.println("PARENT of 3: "+ sets.find_set(3));
//                            System.out.println("PARENT of 4: "+ sets.find_set(4));
//
//                            System.out.println(AdjList[0]);
//                            System.out.println(AdjList[1]);
//                            System.out.println(AdjList[2]);
//                            System.out.println(AdjList[3]);
//                            System.out.println(AdjList[4]);

                            // all src connected
                            unconnected = false; // break;
                        }
                    }
                    // if not unique, then it's duplicated random edge --> continue to next loop
                    else continue;
                }

                // Generating Random Connected Graph COMPLETED ---------------------------------------------------------

                /**
                 * Below gives the initialization of graph_c as a object of MST with size as pickedEdges, the ArrayList
                 * storing all edges.
                 */
                int E = pickedEdges.size();
                MST graph_c = new MST(N, E);
//                System.out.println("here we have new MST with E: "+ E);
                for(int e = 0; e < E; e++){
                    graph_c.edges[e].src = pickedEdges.get(e).src;
                    graph_c.edges[e].dest = pickedEdges.get(e).dest;
                    graph_c.edges[e].weight = pickedEdges.get(e).weight;
                }

                // Timing for Q1_c (Kruskal's algorithm)
                long K_startTime_Q1c = System.nanoTime();
                K_testWeight_Q1c += graph_c.KruskalMST();
                long K_endTime_Q1c = System.nanoTime();
                K_testTime_Q1c += (K_endTime_Q1c - K_startTime_Q1c);

                // Timing for Q1_c (Prim's algorithm)
                long P_startTime_Q1c = System.nanoTime();
                P_testWeight_Q1c += graph_c.PrimMST();
                long P_endTime_Q1c = System.nanoTime();
                P_testTime_Q1c += (P_endTime_Q1c - P_startTime_Q1c);
            }

            K_averWeight_Q1c[s] = K_testWeight_Q1c/50;
            K_averRunTime_Q1c[s] =  K_testTime_Q1c/50;

            P_averWeight_Q1c[s] = P_testWeight_Q1c/50;
            P_averRunTime_Q1c[s] =  P_testTime_Q1c/50;

            System.out.println("Random Connected Graph Size: "+ sizes[s] + " Kruskal's average RunTime: "+ K_averRunTime_Q1c[s]);
            System.out.println("Random Connected Graph Size: "+ sizes[s] + " Kruskal's average Weight: "+ K_averWeight_Q1c[s]);

            System.out.println("Random Connected Graph Size: "+ sizes[s] + " Prim's average RunTime: "+ P_averRunTime_Q1c[s]);
            System.out.println("Random Connected Graph Size: "+ sizes[s] + " Prim's average Weight: "+ P_averWeight_Q1c[s]);
        }


//        for(int s =0; s<sizes.length;s++){
//            System.out.println("Random Connected Graph Size: "+ sizes[s] + "Kruskal's average RunTime: "+ K_averRunTime_Q1c[s]);
//        }
////        for(int s =0; s<sizes.length;s++){
////            System.out.println("Random Connected Graph Size: "+ sizes[s] + "Prim's average RunTime: "+ P_averRunTime_Q1c[s]);
////        }
//        for(int s =0; s<sizes.length;s++){
//            System.out.println("Random Connected Graph Size: "+ sizes[s] + " Kruskal's average Weight: "+ K_averWeight_Q1c[s]);
//        }
////        for(int s =0; s<sizes.length;s++){
////            System.out.println("Random Connected Graph Size: "+ sizes[s] + " Prim's average Weight: "+ P_averWeight_Q1c[s]);
////        }

    }
}
