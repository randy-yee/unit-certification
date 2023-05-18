read("ExperimentFunctions.py");
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


print("Experiment file 1-1");
{
sigstring = "1-1";

b_ranges =
[
    1,12,1;          \\# best 5,  max 50, likely range []
    1,12,1;          \\# best 4 max 50, likely range []
    1,12,1;          \\# best 6 3 max 50, likely range []
    29,50,2;          \\# best 37 max 400, likely range []
    20,40,2;           \\#best 26 max 400, likely range []
    10,35,2;           \\#best 16 6 max 400, likely range []
    30,45,2;           \\#best 37 max 250,
    75,100,2;          \\#best 88 max 500, likely range [60,160]
    85,120,3;         \\#best 97  max 500, likely range [40,120]
    275, 335, 5;       \\#best 305 max 600 , likely range [300,500]
    130, 175, 4;       \\#best 150 max 600 , likely range [100,300]
    175, 230, 5;       \\#best 205 12 max600 , likely range [100,400]
    450, 700, 20;      \\#best 530 max 3300, likely range [600,800]
    600, 800, 20;   \\# best 700 max1500 likely range [500,900]
    400, 700, 20;   \\#15 best 580 max 1500
    150,5000,150;      \\#ceil(i/3)== 6
    150,5000,150;      \\#ceil(i/3)== 6
    150,5000,150       \\#ceil(i/3)== 6
];

/* old ranges

b_ranges =
[
    1,10,1;          \\#max 100, likely range [1,10]
    1,10,1;          \\#max 100, likely range [1,10,]
    1,10,1;          \\#3 max 100, likely range [1,10]
    30,60,3;          \\#max 400, likely range [30,70]
    5,35,3;           \\#max 200, likely range [5,35]
    4,25,2;           \\#6 max 200, likely range [1,30]
    15,60,3;           \\#max 400, likely range [60,160]
    60,160,8;          \\#max 800, likely range [60,160]
    40,120,10;         \\#9 max 400, likely range [40,120]
    300, 500, 15;       \\#max1500 , likely range [300,500]
    100, 300, 25;       \\#max600 , likely range [100,300]
    100, 400, 30;       \\#12 max600 , likely range [100,400]
    300, 700, 20;      \\#max 3300, likely range [600,800]
    500, 900, 50;     \\#max2600 likely range [500,900]
    350, 700, 50;     \\#15 ceil(i/3)== 5
    150,10000,20;      \\#ceil(i/3)== 6
    150,10000,20;      \\#ceil(i/3)== 6
    150,10000,20       \\#ceil(i/3)== 6
];

*/
}
