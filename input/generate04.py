


read("polygenerator.gp");
default(realbitprecision, 1000);
default(parisizemax, 12G);
{


print("hello");

OUTPUT_FILE_PREFIX = "experiment-polynomials-";
STARTING_DISC_RANGE = 32;
ENDING_DISC_RANGE = 129;
DISC_GAP_SIZE = 16;
DISC_SAMPLE_SIZE = 5;
REAL_EMBEDDINGS = 0;
COMPLEX_EMBEDDINGS = 4;

generate_polynomial_list_bitsize(OUTPUT_FILE_PREFIX, \
                                STARTING_DISC_RANGE, ENDING_DISC_RANGE, \
                                DISC_GAP_SIZE, DISC_SAMPLE_SIZE, \
                                REAL_EMBEDDINGS, COMPLEX_EMBEDDINGS);



}
