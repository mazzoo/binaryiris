#include <sys/types.h>
#include <inttypes.h>
#include <sys/stat.h>
#include <limits.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include <errno.h>
#include <math.h>


typedef struct stats_s
{
  uint32_t min;
  uint32_t max;
  double   mean;
  uint32_t meancount;
  uint32_t clipcount;
} stats_t;

#define BUF_SZ 4096
char buf[BUF_SZ];

void usage(void)
{
  printf("usage:\n");
  printf("\tbinaryiris image > iris_of_image.pbm\n");
  printf("\n");
  printf("\tdumps a characteristic 256x256 image of image\n");
  printf("\tthe output is created by sweeping over image and collecting bytewise neighbour relationships\n");
  printf("\n");
  printf("\te.g.:\n");
  printf("\t\t$ binaryiris image > iris_of_image.pbm\n");
  printf("\n");
}

void usage_exit(char * err)
{
  printf("\nERROR: ");
  printf(err);
  printf("\n\n");
  usage();
  exit(1);
}

void gather_stats(uint32_t iris[256][256][3], int color, stats_t * s)
{
  int x;
  int y;

  s->min       = 0xffffffff;
  s->max       = 0;
  s->mean      = 0.0;
  s->meancount = 0;
  s->clipcount = 0;

  for (y = 0; y < 256; y++)
  {
    for (x = 0; x < 256; x++)
    {
      if (iris[x][y][color])
      {
        s->mean += iris[x][y][color];
        s->meancount++;
      }
      if (iris[x][y][color] > s->max)
      {
        s->max = iris[x][y][color];
      }
      if (iris[x][y][color] < s->min)
      {
        s->min = iris[x][y][color];
      }
      if (iris[x][y][color] > 255)
      {
        s->clipcount++;
      }
    }
  }
  s->mean /= s->meancount;
}

void print_stats(stats_t * s)
{
  printf("##############\n");
  printf("# min       : %u\n", s->min);
  printf("# max       : %u\n", s->max);
  printf("# mean      : %lf\n", s->mean);
  printf("# clipcount : %u\n", s->clipcount);
}

int main(int argc, char ** argv)
{
  int ret = 0;

  if (argc != 2)
  {
    usage_exit("one arguments required");
  }

  /* read in file */
  int fimage;
  int fsize;
  uint8_t * image;

  fimage = open(argv[1], O_RDONLY);
  if (fimage < 0)
  {
    snprintf(buf, BUF_SZ, "cannot open file %s", argv[1]);
    usage_exit(buf);
  }
  fsize = lseek(fimage, 0, SEEK_END);
  if (fimage < 0)
  {
    snprintf(buf, BUF_SZ, "cannot seek file %s", argv[1]);
    usage_exit(buf);
  }
  ret = lseek(fimage, 0, SEEK_SET);
  if (ret < 0)
  {
    snprintf(buf, BUF_SZ, "cannot rewind file %s", argv[1]);
    usage_exit(buf);
  }

  image = malloc(fsize);
  if (!image)
  {
    usage_exit("not enough memory");
  }

  ret = read(fimage, image, fsize);
  if (ret != fsize)
  {
    snprintf(buf, BUF_SZ, "read: expected %d but got %d", fsize, ret);
    usage_exit(buf);
  }
  close(fimage);

  int x, y;

  /* PBM header */
  printf("P3\n");
  printf("256 256\n");
  printf("255\n");

  uint32_t iris[256][256][3];
  memset(iris, 0, 3 * 256 * 256 * sizeof(uint32_t));

  /* gather */
  int    i;
  stats_t s0[3];
  stats_t s1[3];
  stats_t s2[3];
  stats_t s3[3];

  /* RED = byte wise neighbours */
  for (i = 1; i < fsize; i++)
  {
    /* suppress 0x00 and 0xff sequences */
    if (
          ( ! ((image[i-1] == 0x00) && (image[i] == 0x00)) ) &&
          ( ! ((image[i-1] == 0xff) && (image[i] == 0xff)) )
       )
    {
      ++iris[image[i-1]][image[i]][0];
    }
  }

  /* GREEN = 16bit wise neighbours */
  for (i = 2; i < fsize; i++)
  {
    /* suppress 0x00 and 0xff sequences */
    if (
          ( ! ((image[i-2] == 0x00) && (image[i] == 0x00)) ) &&
          ( ! ((image[i-2] == 0xff) && (image[i] == 0xff)) )
       )
    {
      ++iris[image[i-2]][image[i]][1];
    }
  }

  /* BLUE = 32bit wise neighbours */
  for (i = 4; i < fsize; i++)
  {
    /* suppress 0x00 and 0xff sequences */
    if (
          ( ! ((image[i-4] == 0x00) && (image[i] == 0x00)) ) &&
          ( ! ((image[i-4] == 0xff) && (image[i] == 0xff)) )
       )
    {
      ++iris[image[i-4]][image[i]][2];
    }
  }

  for (i = 0; i < 3; i++)
  {
    gather_stats(iris, i, &s0[i]);
    print_stats(&s0[i]);
  }

  /* eliminate a minval and shift all to black */
  for (i = 0; i < 3; i++)
  {
    if (s0[i].min)
    {
      for (y = 0; y < 256; y++)
      {
        for (x = 0; x < 256; x++)
        {
          iris[x][y][i] -= s0[i].min;
        }
      }
    }
  }

  for (i = 0; i < 3; i++)
  {
    gather_stats(iris, i, &s1[i]);
    print_stats(&s1[i]);
  }

  /* normalize */
  for (i = 0; i < 3; i++)
  {
    for (y = 0; y < 256; y++)
    {
      for (x = 0; x < 256; x++)
      {
        iris[x][y][i] *= (double) 255.0 / (double) s1[i].max;
      }
    }
  }

  for (i = 0; i < 3; i++)
  {
    gather_stats(iris, i, &s2[i]);
    print_stats(&s2[i]);
  }

#if 1
  double hyperbol[256][3];

  for (i = 0; i < 3; i++)
  {
    s2[i].mean *= s2[i].meancount/(255.0*127.0);

    int j;
    for (j = 0; j < 256; j++)
    {
      if (j < s2[i].mean)
      {
        hyperbol[j][i] = j * 255.0 / s2[i].mean;
      }else{
        hyperbol[j][i] = j * s2[i].mean / 255.0 + 255.0 - s2[i].mean;
      }
    }

    for (y = 0; y < 255; y++)
    {
      for (x = 0; x < 255; x++)
      {
        iris[x][y][i] = hyperbol[iris[x][y][i]][i];
      }
    }
  }
#endif

  for (i = 0; i < 3; i++)
  {
    gather_stats(iris, i, &s3[i]);
    print_stats(&s3[i]);
  }

  /* dump */
  for (y = 0; y < 256; y++)
  {
    for (x = 0; x < 256; x++)
    {
      for (i = 0; i < 3; i++)
      {
        printf("%d ", iris[x][y][i]);
      }
    }
    printf("\n");
  }

  return 0;
}
