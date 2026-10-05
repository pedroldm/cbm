/*
 * Run an external solver (LKH / Linkern) without a shell.
 *
 * Replaces system(): no quoting issues, a real exit status, and the child stays
 * in the caller's process group, so the experiment runner's group kill on a
 * hard timeout also reaches it. Thread-safe (posix_spawn), so CBMLKH's OpenMP
 * workers may call it concurrently.
 */
#ifndef CBM_SPAWN_H
#define CBM_SPAWN_H

#include <errno.h>
#include <fcntl.h>
#include <spawn.h>
#include <sys/types.h>
#include <sys/wait.h>

extern char **environ;

/* Runs argv (argv[0] is the executable path) with stdout and stderr appended to
 * log_path (/dev/null when NULL). Returns the exit status, 128 + signal number
 * if the child was killed, or -1 if it could not be started (errno is set). */
static inline int cbm_run(char *const argv[], const char *log_path)
{
    posix_spawn_file_actions_t actions;
    pid_t pid;
    int rc, status;

    posix_spawn_file_actions_init(&actions);
    posix_spawn_file_actions_addopen(&actions, 0, "/dev/null", O_RDONLY, 0);
    posix_spawn_file_actions_addopen(&actions, 1, log_path ? log_path : "/dev/null", O_WRONLY | O_CREAT | O_APPEND, 0644);
    posix_spawn_file_actions_adddup2(&actions, 1, 2);
    rc = posix_spawn(&pid, argv[0], &actions, NULL, argv, environ);
    posix_spawn_file_actions_destroy(&actions);
    if (rc != 0) {
        errno = rc;
        return -1;
    }

    while (waitpid(pid, &status, 0) < 0) {
        if (errno != EINTR) return -1;
    }
    if (WIFEXITED(status)) return WEXITSTATUS(status);
    if (WIFSIGNALED(status)) return 128 + WTERMSIG(status);
    return -1;
}

#endif /* CBM_SPAWN_H */
