
#include "GA_Func.h"
// Number of individuals in each generation 
#define POPULATION_SIZE 100 //Make this even, I divide by an even number below and it throws out the 
int MaxGenerations = 100;
const int len = 3;
std::array<double,len>  GENES = {10,10,10};
std::array<double,len>  LowerGENES = {0,0,0};
std::array<double,len>  UpperGENES = {200,200,200};


std::default_random_engine generator;


// Create random genes for mutation 
double mutated_genes(int i, double mean, double std) 
{ 
     std::normal_distribution<double> distribution(mean, std);
     double r = distribution(generator);
     if (r<LowerGENES[i]){
         r = LowerGENES[i];
     }
     else if(r > UpperGENES[i]){
         r = UpperGENES[i];
     } 
    return r;
} 


// Perform mating and produce new offspring 
Individual mate(Individual par1, Individual par2) 
{ 
    // chromosome for offspring 
    Individual child_chromosome; 
    double mean;
    double std;
    for(int i = 0;i<len;i++) 
    { 
        // random probability  
        float p = random_num(0, 100)/100; 
        // if prob is less than 0.45, insert gene 
        // from parent 1  
        if(p < 0.45) {
            child_chromosome.chromosome.push_back(par1.chromosome[i]); 
            //std::cout << "Parent 1" << std::endl;
        }
        // if prob is between 0.45 and 0.90, insert 
        // gene from parent 2 
        else if(p < 0.90){
            child_chromosome.chromosome.push_back(par2.chromosome[i]); 
            //std::cout << "Parent 2" << std::endl;
        }
        // otherwise insert random gene(mutate),  
        // for maintaining diversity 
        else{
            //std::cout << "Mutation" << std::endl;
            mean = (par1.chromosome[i] + par2.chromosome[i]) /2;
            std = (UpperGENES[i]-LowerGENES[i])*0.1; // One standard deviation is 10% of the size of upper-lower bound
            child_chromosome.chromosome.push_back(mutated_genes(i,mean,std)); 
        }
    } 
    //child_chromosome.fitness = fitness_function(child_chromosome.chromosome);
    //std::cout << "Parent 1: "<< par1.chromosome[0] << " " <<par1.chromosome[1] << " " <<par1.chromosome[2] << std::endl;
    //std::cout << "Parent 2: "<< par2.chromosome[0] << " " <<par2.chromosome[1] << " " <<par2.chromosome[2] << std::endl;
    //std::cout << "child_chromosome: "<< child_chromosome.chromosome[0] << " " <<child_chromosome.chromosome[1] << " " <<child_chromosome.chromosome[2] << std::endl;
    
    return child_chromosome; 
}; 

std::vector<double> create_gnome() 
{ 
    std::vector<double> gnome;
    for(int i = 0;i<len;i++){
        gnome.push_back(random_num(LowerGENES[i], UpperGENES[i]));
        //std::cout << "LOOP" << std::endl;
    }
    return gnome; 
} 


bool CustomSortFitness(Individual const& lhs, Individual const& rhs) {
        return lhs.fitness < rhs.fitness;
}

//#########################################     Multi Threading Implementation          #########################################
//###############################################################################################################################
std::mutex mtx_pop,mtx_newgen;
#define NUM_THREADS 4 //thread::hardware_concurrency();
std::stack<Individual> WorkerStack;
std::condition_variable cv;
std::condition_variable cv_main;
bool Processed = false;


void WorkerPool() {
   Individual result1;

   while(1){ 
       
        std::unique_lock<std::mutex> lck(mtx_pop);
        //std::cout << "Waiting for Command\n";
        cv.wait(lck,[]{return !WorkerStack.empty() && Processed;});
        //std::cout <<  "Worker Stack Size: " << WorkerStack.size() << "\n" << std::flush;
        result1 = WorkerStack.top();
        WorkerStack.pop();
        //std::cout << "InThreads\t" << WorkerStack.size() << "\n";
        //std::cout << result1.chromosome[0] <<  "\t" << result1.chromosome [1] << "\t" << result1.chromosome [1] << "\n";
        lck.unlock();
        
        
    //Calculate the fitness of the individual you pulled off
    result1.fitness = fitness_function(result1.chromosome);
    //Mutex to put the Individiual back into the new_generation after calculating the fitness function
        //std::cout << "Trying to get lock\n";
        std::unique_lock<std::mutex> lck_newgen(mtx_newgen);
        //std::cout << "GOT LOCK\n";
        //std::cout <<  "Thread: " << "\t" << result1.chromosome[0] <<  "\t" << result1.chromosome [1] << "\t" << result1.chromosome [1] << "\t" << result1.fitness << "\n";
        //std::cout << "Done Threads\n";
        new_generation.push_back(result1);
        //std::cout <<  "New Gen Size: " << new_generation.size() << "\n" << std::flush;
        lck_newgen.unlock();
        if (new_generation.size() == POPULATION_SIZE){
            cv_main.notify_one();
        }
   }
}



int main(int argc, char* argv[]) { 
    int SatsPerPlane = 5;
    int Planes = 2;
    satellite SeedSat = CreateSat ( 14.7335, 0, 6986, 0, 0, 186.606);
    constellation WalkerConst = CreateConstellation  (SeedSat, SatsPerPlane, Planes);
    exit(-1);



   std::cout << "argc = " << argc << std::endl; 
   for(int i = 0; i < argc; i++) 
      std::cout << "argv[" << i << "] = " << argv[i] << std::endl;

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();   

    srand((unsigned)(time(0)));
    //Creating thread pool
    int rc,i;
    std::vector<std::thread> Pool;
    for( i = 0; i < NUM_THREADS; i++ ) {
      //std::cout << "main() : creating thread, " << i << std::endl;
      Pool.push_back(std::thread(WorkerPool));
      if (rc) {
         std::cout << "Error:unable to create thread," << rc << std::endl;
         exit(-1);
      }
   }
    
    // current generation 
    int generation = 0; 
  
    for(int i = 0;i<POPULATION_SIZE;i++) 
    { 
        Individual gnome;
        gnome.chromosome = create_gnome(); 
        WorkerStack.push(gnome);
    } 
    
    //Wait for threads to finish up calculating all the fitness values and populating the new_generation vector
    //std::cout << "OG Size: " << WorkerStack.size() << "\n";
    //std::cout << "ABOUT TO unlock main mutex";
    //std::cout << "setting processed to true"  << std::flush;;
    Processed = true;
    cv.notify_all();
    //std::cout << "Unlocked main mutex"  << std::flush;;
    std::unique_lock<std::mutex> lck_mnewgen(mtx_newgen);
    //std::cout << "GEN SIZE: " << new_generation.size();
    cv_main.wait(lck_mnewgen,[]{return new_generation.size() == POPULATION_SIZE;});
    //std::cout << "DONE\n"  << std::flush;;
    population = new_generation;
    Processed = false;
    
 //   exit(-1);
    while(generation < MaxGenerations){
        //sort population by fitness value
        //std::cout << population.begin(); 
        //Start creating the new generation
        std::vector<Individual> ().swap(new_generation);
        //Perform elietism, top 10% move onto the next generation
        //int s = (10*POPULATION_SIZE)/100;;
        int s = ((10*POPULATION_SIZE) / 100) + (((10*POPULATION_SIZE) % 100) != 0);
        //std::cout <<"Ceiling: " << s;
        for(int i = 0;i<s;i++)
            new_generation.push_back(population[i]); 

        // for the next 90% have the parents mate to create the next generation
        s = (90*POPULATION_SIZE)/100;
        for(int i = 0;i<s;i++) 
        { 
            int r = random_num(0, POPULATION_SIZE/2); //Only using the best half of the population as the future parents
            Individual parent1 = population[r]; 
            r = random_num(0, POPULATION_SIZE/2);     //Only using the best half of the population as the future parents
            Individual parent2 = population[r]; 
            Individual offspring = mate(parent1, parent2); 
            WorkerStack.push(offspring);
            //std::cout << "Worker Stack Main: " << WorkerStack.size();
            //new_generation.push_back(offspring);  
        }
        
        Processed = true;
        cv.notify_all();
        //std::cout << "\nDone With Gen" << generation << "\n";
        cv_main.wait(lck_mnewgen,[]{return new_generation.size() == POPULATION_SIZE;});
        //std::cout <<"Done Waiting for new gen \n";
        population = new_generation;
        Processed = false;
        sort(population.begin(), population.end(), &CustomSortFitness) ;

        std::cout<< "Generation: " << generation << "\t"; 
        std::cout<< "Array Values: "<< population[0].chromosome[0] << " " << population[0].chromosome[1] << " " << population[0].chromosome[2] <<"\t"; 
        std::cout<< "Fitness: "<< population[0].fitness << "\n"; 
        generation++;
    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);
    std::cout << "It took me " << time_span.count() << " seconds.";
    std::cout << std::endl;  
}