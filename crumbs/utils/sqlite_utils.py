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

    def __getitem__(self, key):
        c = self.connection.cursor()
        sql = "select * from {cache_table} where key = '{key}'"
        c.execute(sql.format(key=key, cache_table=self._cache_table_name))

        result = c.fetchall()
        if not result:
            return None
        if len(result) > 1:
            msg = 'More than one ocurrence of the key:{}'.format(key)
            raise RuntimeError(msg)

        return pickle.loads(str(result[0][1]))

    def __setitem__(self, key, value):
        c = self.connection.cursor()
        exists = self[key]
        value = pickle.dumps(value, pickle.HIGHEST_PROTOCOL)
        if exists:
            sql = "update {cache_table} set value = ? where key = ?"
            sql = sql.format(cache_table=self._cache_table_name, key=key)
            c.execute(sql, [sqlite3.Binary(value), key])
        else:
            sql = "insert into {cache_table} values(?, ?)"
            sql = sql.format(cache_table=self._cache_table_name)
            c.execute(sql, [key, sqlite3.Binary(value)])
        self.connection.commit()

    def close(self):
        if self._fhand is not None:
            self._fhand.close()
        self.connection.close()

    def __str__(self):
        c = self.connection.cursor()
        text = None
        sql = "select * from {}".format(self._cache_table_name)
        for key, value in c.execute(sql):
            value = pickle.loads(str(value))
            if text is None:
                key_order = value.keys()
                text = ['\t'.join(['key'] + key_order)]
            text.append('\t'.join([key] + [str(value[valkey])
                                           for valkey in key_order]))
        return '\n'.join(text)
