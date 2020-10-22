#include <stdio.h>

int main(void)
{
    int N, num_of_steps;
    scanf("%d%d", &N, &num_of_steps);
    int *line_numbers = new int[num_of_steps];
    
    FILE* fid_lines = fopen("lines.txt", "r");
    for (int n = 0; n < num_of_steps; ++n)
    {
        fscanf(fid_lines, "%d", &line_numbers[n]);
        line_numbers[n] += 2;
    }
    fclose(fid_lines);

    FILE* fid_in = fopen("OUTCAR", "r");
    FILE* fid_out = fopen("xf.txt", "w");
    char line[256];
    int i = 0;
    int n = 0;
    while (fgets(line, sizeof(line), fid_in))
    {
        i++;
        if(i >= line_numbers[n] && i < line_numbers[n] + N)
        {
            fprintf(fid_out, "%s", line);
            if (i == line_numbers[n] + N - 1) { n++; }
        }
    }
    fclose(fid_in);
    fclose(fid_out);
    delete[] line_numbers;
    return 0;
}

