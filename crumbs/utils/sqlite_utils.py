'''
Created on 2014 aza 3

@author: peio
'''
import os
import sqlite3
import cPickle as pickle


class SqliteCache(object):
    def __init__(self, fpath):
        self._cache_table_name = 'cachedata'
        self._fhand = None
        if not os.path.exists(fpath):
            self._fhand = open(fpath, 'w')
            fpath = self._fhand.name

        self.connection = sqlite3.connect(fpath)
        self._create_init_table()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()
        if type:
            return 0

    def _create_init_table(self):
        c = self.connection.cursor()
        c.execute("SELECT name FROM sqlite_master WHERE type='table';")
        result = c.fetchall()
        if not result or self._cache_table_name not in result[0]:
            sql = 'create table {table_name} (key TEXT, value BLOB)'
            c.execute(sql.format(table_name=self._cache_table_name))
        self.connection.commit()

    def get(self, key):
        c = self.connection.cursor()
        sql = "select * from {cache_table} where key = '{key}'"
        c.execute(sql.format(key=key, cache_table=self._cache_table_name))
        result = c.fetchone()

        if result is not None:
            result = pickle.loads(str(result[1]))
        return result

    def set(self, key, value):
        c = self.connection.cursor()
        value = pickle.dumps(value, pickle.HIGHEST_PROTOCOL)
        sql = "insert into {cache_table} values(?, ?)"
        sql = sql.format(cache_table=self._cache_table_name)
        c.execute(sql, [key, sqlite3.Binary(value)])
        self.connection.commit()

    def close(self):
        if self._fhand is not None:
            self._fhand.close()
        self.connection.close()
