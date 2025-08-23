


read("polygenerator.gp");
default(realbitprecision, 1000);
default(parisizemax, 12G);
{


print("hello");

OUTPUT_FILE_PREFIX = "new-experiment-polynomials-";
STARTING_REG_RANGE = 16; \\ ~32 disc
ENDING_REG_RANGE = 65;   \\ ~128-ish
REG_GAP_SIZE = 8;        \\ ~disc gap 16
SAMPLE_SIZE = 5;
REAL_EMBEDDINGS = 4;
COMPLEX_EMBEDDINGS = 0;

generate_polynomial_list_bitsize(OUTPUT_FILE_PREFIX, \
                                STARTING_REG_RANGE, ENDING_REG_RANGE, \
                                REG_GAP_SIZE, SAMPLE_SIZE, \
                                REAL_EMBEDDINGS, COMPLEX_EMBEDDINGS);



}
