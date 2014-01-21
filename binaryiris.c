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

void gather_stats(uint32_t iris[256][256], stats_t * s)
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
      if (iris[x][y])
      {
        s->mean += iris[x][y];
        s->meancount++;
      }
      if (iris[x][y] > s->max)
      {
        s->max = iris[x][y];
      }
      if (iris[x][y] < s->min)
      {
        s->min = iris[x][y];
      }
      if (iris[x][y] > 255)
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
  printf("P2\n");
  printf("256 256\n");
  printf("256\n");

  uint32_t iris[256][256];
  memset(iris, 0, 256 * 256 * sizeof(uint32_t));

  /* gather */
  int    i;
  stats_t s0;
  stats_t s1;
  stats_t s2;
  stats_t s3;

  for (i = 1; i < fsize; i++)
  {
    /* suppress 0x00 and 0xff sequences */
    if (
          ( ! ((image[i-1] == 0x00) && (image[i] == 0x00)) ) &&
          ( ! ((image[i-1] == 0xff) && (image[i] == 0xff)) )
       )
    {
      ++iris[image[i-1]][image[i]];
    }
  }

  gather_stats(iris, &s0);
  print_stats(&s0);

  /* eliminate a minval and shift all to black */
  if (s0.min)
  {
    for (y = 0; y < 256; y++)
    {
      for (x = 0; x < 256; x++)
      {
        iris[x][y] -= s0.min;
      }
    }
  }

  gather_stats(iris, &s1);
  print_stats(&s1);

  /* normalize */
  for (y = 0; y < 256; y++)
  {
    for (x = 0; x < 256; x++)
    {
      iris[x][y] *= (double) 255.0 / (double) s1.max;
    }
  }

  gather_stats(iris, &s2);
  print_stats(&s2);

  double hyperbol[256];

  s2.mean *= 2;
  for (i = 0; i < 256; i++)
  {
    if (i < s2.mean)
    {
      hyperbol[i] = i * 255.0 / s2.mean;
    }else{
      hyperbol[i] = i * s2.mean / 255.0 + 255.0 - s2.mean;
    }
  }

  for (y = 0; y < 255; y++)
  {
    for (x = 0; x < 255; x++)
    {
      iris[x][y] = hyperbol[iris[x][y]];
    }
  }

  gather_stats(iris, &s3);
  print_stats(&s3);

  /* dump */
  for (y = 0; y < 256; y++)
  {
    for (x = 0; x < 256; x++)
    {
      printf("%d ", iris[x][y]);
    }
    printf("\n");
  }

  return 0;
}
