/*
  File autogenerated by gengetopt version 2.22.4
  generated with the following command:
  /genome4/raid/gus/Source/gengetopt-2.22.4/src/gengetopt --show-required --default-option -c cc -H hpp -F hash_merge_cmdline -f hash_merge_cmdline -a hash_merge_args --unamed-opts=database.jf

  The developers of gengetopt consider the fixed text that goes in all
  gengetopt output files to be in the public domain:
  we make no copyright claims on it.
*/

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef FIX_UNUSED
#define FIX_UNUSED(X) (void) (X) /* avoid warnings for unused params */
#endif

#include <getopt.h>

#include "hash_merge_cmdline.hpp"

const char *hash_merge_args_purpose = "Merge jellyfish databases";

const char *hash_merge_args_usage = "Usage: jellyfish merge [OPTIONS]... [database.jf]...";

const char *hash_merge_args_description = "";

const char *hash_merge_args_help[] = {
  "  -h, --help                    Print help and exit",
  "  -V, --version                 Print version and exit",
  "  -s, --buffer-size=Buffer length\n                                Length in bytes of input buffer  \n                                  (default=`10000000')",
  "  -o, --output=STRING           Output file  (default=`mer_counts_merged.jf')",
  "      --out-counter-len=INT     Length (in bytes) of counting field in output  \n                                  (default=`4')",
  "      --out-buffer-size=LONG    Size of output buffer per thread  \n                                  (default=`10000000')",
  "  -v, --verbose                 Be verbose  (default=off)",
    0
};

typedef enum {ARG_NO
  , ARG_FLAG
  , ARG_STRING
  , ARG_INT
  , ARG_LONG
} hash_merge_cmdline_arg_type;

static
void clear_given (struct hash_merge_args *args_info);
static
void clear_args (struct hash_merge_args *args_info);

static int
hash_merge_cmdline_internal (int argc, char **argv, struct hash_merge_args *args_info,
                        struct hash_merge_cmdline_params *params, const char *additional_error);


static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct hash_merge_args *args_info)
{
  args_info->help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->buffer_size_given = 0 ;
  args_info->output_given = 0 ;
  args_info->out_counter_len_given = 0 ;
  args_info->out_buffer_size_given = 0 ;
  args_info->verbose_given = 0 ;
}

static
void clear_args (struct hash_merge_args *args_info)
{
  FIX_UNUSED (args_info);
  args_info->buffer_size_arg = 10000000;
  args_info->buffer_size_orig = NULL;
  args_info->output_arg = gengetopt_strdup ("mer_counts_merged.jf");
  args_info->output_orig = NULL;
  args_info->out_counter_len_arg = 4;
  args_info->out_counter_len_orig = NULL;
  args_info->out_buffer_size_arg = 10000000;
  args_info->out_buffer_size_orig = NULL;
  args_info->verbose_flag = 0;
  
}

static
void init_args_info(struct hash_merge_args *args_info)
{


  args_info->help_help = hash_merge_args_help[0] ;
  args_info->version_help = hash_merge_args_help[1] ;
  args_info->buffer_size_help = hash_merge_args_help[2] ;
  args_info->output_help = hash_merge_args_help[3] ;
  args_info->out_counter_len_help = hash_merge_args_help[4] ;
  args_info->out_buffer_size_help = hash_merge_args_help[5] ;
  args_info->verbose_help = hash_merge_args_help[6] ;
  
}

void
hash_merge_cmdline_print_version (void)
{
  printf ("%s %s\n",
     (strlen(HASH_MERGE_CMDLINE_PACKAGE_NAME) ? HASH_MERGE_CMDLINE_PACKAGE_NAME : HASH_MERGE_CMDLINE_PACKAGE),
     HASH_MERGE_CMDLINE_VERSION);
}

static void print_help_common(void) {
  hash_merge_cmdline_print_version ();

  if (strlen(hash_merge_args_purpose) > 0)
    printf("\n%s\n", hash_merge_args_purpose);

  if (strlen(hash_merge_args_usage) > 0)
    printf("\n%s\n", hash_merge_args_usage);

  printf("\n");

  if (strlen(hash_merge_args_description) > 0)
    printf("%s\n\n", hash_merge_args_description);
}

void
hash_merge_cmdline_print_help (void)
{
  int i = 0;
  print_help_common();
  while (hash_merge_args_help[i])
    printf("%s\n", hash_merge_args_help[i++]);
}

void
hash_merge_cmdline_init (struct hash_merge_args *args_info)
{
  clear_given (args_info);
  clear_args (args_info);
  init_args_info (args_info);

  args_info->inputs = 0;
  args_info->inputs_num = 0;
}

void
hash_merge_cmdline_params_init(struct hash_merge_cmdline_params *params)
{
  if (params)
    { 
      params->override = 0;
      params->initialize = 1;
      params->check_required = 1;
      params->check_ambiguity = 0;
      params->print_errors = 1;
    }
}

struct hash_merge_cmdline_params *
hash_merge_cmdline_params_create(void)
{
  struct hash_merge_cmdline_params *params = 
    (struct hash_merge_cmdline_params *)malloc(sizeof(struct hash_merge_cmdline_params));
  hash_merge_cmdline_params_init(params);  
  return params;
}

static void
free_string_field (char **s)
{
  if (*s)
    {
      free (*s);
      *s = 0;
    }
}


static void
hash_merge_cmdline_release (struct hash_merge_args *args_info)
{
  unsigned int i;
  free_string_field (&(args_info->buffer_size_orig));
  free_string_field (&(args_info->output_arg));
  free_string_field (&(args_info->output_orig));
  free_string_field (&(args_info->out_counter_len_orig));
  free_string_field (&(args_info->out_buffer_size_orig));
  
  
  for (i = 0; i < args_info->inputs_num; ++i)
    free (args_info->inputs [i]);

  if (args_info->inputs_num)
    free (args_info->inputs);

  clear_given (args_info);
}


static void
write_into_file(FILE *outfile, const char *opt, const char *arg, const char *values[])
{
  FIX_UNUSED (values);
  if (arg) {
    fprintf(outfile, "%s=\"%s\"\n", opt, arg);
  } else {
    fprintf(outfile, "%s\n", opt);
  }
}


int
hash_merge_cmdline_dump(FILE *outfile, struct hash_merge_args *args_info)
{
  int i = 0;

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot dump options to stream\n", HASH_MERGE_CMDLINE_PACKAGE);
      return EXIT_FAILURE;
    }

  if (args_info->help_given)
    write_into_file(outfile, "help", 0, 0 );
  if (args_info->version_given)
    write_into_file(outfile, "version", 0, 0 );
  if (args_info->buffer_size_given)
    write_into_file(outfile, "buffer-size", args_info->buffer_size_orig, 0);
  if (args_info->output_given)
    write_into_file(outfile, "output", args_info->output_orig, 0);
  if (args_info->out_counter_len_given)
    write_into_file(outfile, "out-counter-len", args_info->out_counter_len_orig, 0);
  if (args_info->out_buffer_size_given)
    write_into_file(outfile, "out-buffer-size", args_info->out_buffer_size_orig, 0);
  if (args_info->verbose_given)
    write_into_file(outfile, "verbose", 0, 0 );
  

  i = EXIT_SUCCESS;
  return i;
}

int
hash_merge_cmdline_file_save(const char *filename, struct hash_merge_args *args_info)
{
  FILE *outfile;
  int i = 0;

  outfile = fopen(filename, "w");

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot open file for writing: %s\n", HASH_MERGE_CMDLINE_PACKAGE, filename);
      return EXIT_FAILURE;
    }

  i = hash_merge_cmdline_dump(outfile, args_info);
  fclose (outfile);

  return i;
}

void
hash_merge_cmdline_free (struct hash_merge_args *args_info)
{
  hash_merge_cmdline_release (args_info);
}

/** @brief replacement of strdup, which is not standard */
char *
gengetopt_strdup (const char *s)
{
  char *result = 0;
  if (!s)
    return result;

  result = (char*)malloc(strlen(s) + 1);
  if (result == (char*)0)
    return (char*)0;
  strcpy(result, s);
  return result;
}

int
hash_merge_cmdline (int argc, char **argv, struct hash_merge_args *args_info)
{
  return hash_merge_cmdline2 (argc, argv, args_info, 0, 1, 1);
}

int
hash_merge_cmdline_ext (int argc, char **argv, struct hash_merge_args *args_info,
                   struct hash_merge_cmdline_params *params)
{
  int result;
  result = hash_merge_cmdline_internal (argc, argv, args_info, params, 0);

  if (result == EXIT_FAILURE)
    {
      hash_merge_cmdline_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
hash_merge_cmdline2 (int argc, char **argv, struct hash_merge_args *args_info, int override, int initialize, int check_required)
{
  int result;
  struct hash_merge_cmdline_params params;
  
  params.override = override;
  params.initialize = initialize;
  params.check_required = check_required;
  params.check_ambiguity = 0;
  params.print_errors = 1;

  result = hash_merge_cmdline_internal (argc, argv, args_info, &params, 0);

  if (result == EXIT_FAILURE)
    {
      hash_merge_cmdline_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
hash_merge_cmdline_required (struct hash_merge_args *args_info, const char *prog_name)
{
  FIX_UNUSED (args_info);
  FIX_UNUSED (prog_name);
  return EXIT_SUCCESS;
}


static char *package_name = 0;

/**
 * @brief updates an option
 * @param field the generic pointer to the field to update
 * @param orig_field the pointer to the orig field
 * @param field_given the pointer to the number of occurrence of this option
 * @param prev_given the pointer to the number of occurrence already seen
 * @param value the argument for this option (if null no arg was specified)
 * @param possible_values the possible values for this option (if specified)
 * @param default_value the default value (in case the option only accepts fixed values)
 * @param arg_type the type of this option
 * @param check_ambiguity @see hash_merge_cmdline_params.check_ambiguity
 * @param override @see hash_merge_cmdline_params.override
 * @param no_free whether to free a possible previous value
 * @param multiple_option whether this is a multiple option
 * @param long_opt the corresponding long option
 * @param short_opt the corresponding short option (or '-' if none)
 * @param additional_error possible further error specification
 */
static
int update_arg(void *field, char **orig_field,
               unsigned int *field_given, unsigned int *prev_given, 
               char *value, const char *possible_values[],
               const char *default_value,
               hash_merge_cmdline_arg_type arg_type,
               int check_ambiguity, int override,
               int no_free, int multiple_option,
               const char *long_opt, char short_opt,
               const char *additional_error)
{
  char *stop_char = 0;
  const char *val = value;
  int found;
  char **string_field;
  FIX_UNUSED (field);

  stop_char = 0;
  found = 0;

  if (!multiple_option && prev_given && (*prev_given || (check_ambiguity && *field_given)))
    {
      if (short_opt != '-')
        fprintf (stderr, "%s: `--%s' (`-%c') option given more than once%s\n", 
               package_name, long_opt, short_opt,
               (additional_error ? additional_error : ""));
      else
        fprintf (stderr, "%s: `--%s' option given more than once%s\n", 
               package_name, long_opt,
               (additional_error ? additional_error : ""));
      return 1; /* failure */
    }

  FIX_UNUSED (default_value);
    
  if (field_given && *field_given && ! override)
    return 0;
  if (prev_given)
    (*prev_given)++;
  if (field_given)
    (*field_given)++;
  if (possible_values)
    val = possible_values[found];

  switch(arg_type) {
  case ARG_FLAG:
    *((int *)field) = !*((int *)field);
    break;
  case ARG_INT:
    if (val) *((int *)field) = strtol (val, &stop_char, 0);
    break;
  case ARG_LONG:
    if (val) *((long *)field) = (long)strtol (val, &stop_char, 0);
    break;
  case ARG_STRING:
    if (val) {
      string_field = (char **)field;
      if (!no_free && *string_field)
        free (*string_field); /* free previous string */
      *string_field = gengetopt_strdup (val);
    }
    break;
  default:
    break;
  };

  /* check numeric conversion */
  switch(arg_type) {
  case ARG_INT:
  case ARG_LONG:
    if (val && !(stop_char && *stop_char == '\0')) {
      fprintf(stderr, "%s: invalid numeric value: %s\n", package_name, val);
      return 1; /* failure */
    }
    break;
  default:
    ;
  };

  /* store the original value */
  switch(arg_type) {
  case ARG_NO:
  case ARG_FLAG:
    break;
  default:
    if (value && orig_field) {
      if (no_free) {
        *orig_field = value;
      } else {
        if (*orig_field)
          free (*orig_field); /* free previous string */
        *orig_field = gengetopt_strdup (value);
      }
    }
  };

  return 0; /* OK */
}


int
hash_merge_cmdline_internal (
  int argc, char **argv, struct hash_merge_args *args_info,
                        struct hash_merge_cmdline_params *params, const char *additional_error)
{
  int c;	/* Character of the parsed option.  */

  int error = 0;
  struct hash_merge_args local_args_info;
  
  int override;
  int initialize;
  int check_required;
  int check_ambiguity;
  
  package_name = argv[0];
  
  override = params->override;
  initialize = params->initialize;
  check_required = params->check_required;
  check_ambiguity = params->check_ambiguity;

  if (initialize)
    hash_merge_cmdline_init (args_info);

  hash_merge_cmdline_init (&local_args_info);

  optarg = 0;
  optind = 0;
  opterr = params->print_errors;
  optopt = '?';

  while (1)
    {
      int option_index = 0;

      static struct option long_options[] = {
        { "help",	0, NULL, 'h' },
        { "version",	0, NULL, 'V' },
        { "buffer-size",	1, NULL, 's' },
        { "output",	1, NULL, 'o' },
        { "out-counter-len",	1, NULL, 0 },
        { "out-buffer-size",	1, NULL, 0 },
        { "verbose",	0, NULL, 'v' },
        { 0,  0, 0, 0 }
      };

      c = getopt_long (argc, argv, "hVs:o:v", long_options, &option_index);

      if (c == -1) break;	/* Exit from `while (1)' loop.  */

      switch (c)
        {
        case 'h':	/* Print help and exit.  */
          hash_merge_cmdline_print_help ();
          hash_merge_cmdline_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 'V':	/* Print version and exit.  */
          hash_merge_cmdline_print_version ();
          hash_merge_cmdline_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 's':	/* Length in bytes of input buffer.  */
        
        
          if (update_arg( (void *)&(args_info->buffer_size_arg), 
               &(args_info->buffer_size_orig), &(args_info->buffer_size_given),
              &(local_args_info.buffer_size_given), optarg, 0, "10000000", ARG_LONG,
              check_ambiguity, override, 0, 0,
              "buffer-size", 's',
              additional_error))
            goto failure;
        
          break;
        case 'o':	/* Output file.  */
        
        
          if (update_arg( (void *)&(args_info->output_arg), 
               &(args_info->output_orig), &(args_info->output_given),
              &(local_args_info.output_given), optarg, 0, "mer_counts_merged.jf", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "output", 'o',
              additional_error))
            goto failure;
        
          break;
        case 'v':	/* Be verbose.  */
        
        
          if (update_arg((void *)&(args_info->verbose_flag), 0, &(args_info->verbose_given),
              &(local_args_info.verbose_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "verbose", 'v',
              additional_error))
            goto failure;
        
          break;

        case 0:	/* Long option with no short option */
          /* Length (in bytes) of counting field in output.  */
          if (strcmp (long_options[option_index].name, "out-counter-len") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->out_counter_len_arg), 
                 &(args_info->out_counter_len_orig), &(args_info->out_counter_len_given),
                &(local_args_info.out_counter_len_given), optarg, 0, "4", ARG_INT,
                check_ambiguity, override, 0, 0,
                "out-counter-len", '-',
                additional_error))
              goto failure;
          
          }
          /* Size of output buffer per thread.  */
          else if (strcmp (long_options[option_index].name, "out-buffer-size") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->out_buffer_size_arg), 
                 &(args_info->out_buffer_size_orig), &(args_info->out_buffer_size_given),
                &(local_args_info.out_buffer_size_given), optarg, 0, "10000000", ARG_LONG,
                check_ambiguity, override, 0, 0,
                "out-buffer-size", '-',
                additional_error))
              goto failure;
          
          }
          
          break;
        case '?':	/* Invalid option.  */
          /* `getopt_long' already printed an error message.  */
          goto failure;

        default:	/* bug: option not considered.  */
          fprintf (stderr, "%s: option unknown: %c%s\n", HASH_MERGE_CMDLINE_PACKAGE, c, (additional_error ? additional_error : ""));
          abort ();
        } /* switch */
    } /* while */




  hash_merge_cmdline_release (&local_args_info);

  if ( error )
    return (EXIT_FAILURE);

  if (optind < argc)
    {
      int i = 0 ;
      int found_prog_name = 0;
      /* whether program name, i.e., argv[0], is in the remaining args
         (this may happen with some implementations of getopt,
          but surely not with the one included by gengetopt) */

      i = optind;
      while (i < argc)
        if (argv[i++] == argv[0]) {
          found_prog_name = 1;
          break;
        }
      i = 0;

      args_info->inputs_num = argc - optind - found_prog_name;
      args_info->inputs =
        (char **)(malloc ((args_info->inputs_num)*sizeof(char *))) ;
      while (optind < argc)
        if (argv[optind++] != argv[0])
          args_info->inputs[ i++ ] = gengetopt_strdup (argv[optind-1]) ;
    }

  return 0;

failure:
  
  hash_merge_cmdline_release (&local_args_info);
  return (EXIT_FAILURE);
}